#pragma once

#include "lajolla.h"
#include "pcg.h"
#include "phase_function.h"
#include "ray.h"
#include "spectrum.h"
#include "intersection.h"
#include "point_and_normal.h"
#include "light.h"
#include "scene.h"

int update_medium_id(const PathVertex &vertex,
                            const Ray &ray,
                            int medium_id) {
    if (vertex.interior_medium_id != vertex.exterior_medium_id) {
        // At medium transition. Update medium id.
        if (dot(ray.dir, vertex.geometric_normal) > 0) {
            medium_id = vertex.exterior_medium_id;
        } else {
            medium_id = vertex.interior_medium_id;
        }
    }
    return medium_id;
}

Spectrum next_event_estimation(Vector3 pos, Vector3 inc_dir, int med_id, int bounce_cnt, Light src_light, const Scene& scn, pcg32_state& rng_state, int src_id) {

    PointAndNormal light_sample = sample_point_on_light(src_light, pos, 
                                Vector2(next_pcg32_real<Real>(rng_state), next_pcg32_real<Real>(rng_state)),
                                next_pcg32_real<Real>(rng_state), scn);

    Real transmittance = 1;
    int active_medium = med_id;
    int bounces_made = 0;
    Real path_prob = 1;

    Vector3 start_pos = pos;

    while (true) {
        Vector3 ray_dir = normalize(light_sample.position - pos);
        Ray trace_ray = Ray{pos, ray_dir, get_shadow_epsilon(scn),
                            (1 - get_shadow_epsilon(scn)) * distance(light_sample.position, pos)};
        RayDifferential diff_ray = RayDifferential{ Real(0), Real(0) };
        std::optional<PathVertex> hit = intersect(scn, trace_ray, diff_ray);
        Real ray_length = distance(pos, light_sample.position);

        if (hit) {
            ray_length = distance(pos, hit->position);
        }
        if (active_medium != -1) {
            Medium medium_props = scn.media[active_medium];
            Real absorption = get_sigma_a(medium_props, Vector3(1, 2, 3)).x;
            Real scatter = get_sigma_s(medium_props, Vector3(1, 2, 3)).x;
            Real extinction = absorption + scatter;
            Real decay = exp(-extinction * ray_length);
            transmittance *= decay;
            path_prob *= decay;
        }

        if (!hit) break;

        PathVertex vertex = hit.value();
        if (vertex.material_id >= 0) {
            return make_zero_spectrum();
        }

        bounces_made++;
        if (scn.options.max_depth != -1 && bounce_cnt + bounces_made + 1 >= scn.options.max_depth) {
            return make_zero_spectrum();
        }

        active_medium = update_medium_id(vertex, trace_ray, active_medium);
        pos += ray_length * trace_ray.dir;
    }

    if (transmittance > 0) {
        Vector3 outgoing_dir = normalize(start_pos - light_sample.position);
        Real sq_distance = distance_squared(start_pos, light_sample.position);
        Real cos_theta = abs(dot(outgoing_dir, light_sample.normal));
        Real geo_term = cos_theta / sq_distance;

        PhaseFunction phase_func = get_phase_function(scn.media[med_id]);
        Spectrum bsdf = eval(phase_func, inc_dir, -outgoing_dir);
        Spectrum emitted = emission(src_light, outgoing_dir, Real(0), light_sample, scn);
        Real illumination_pdf = light_pmf(scn, src_id) * pdf_point_on_light(src_light, light_sample, start_pos, scn);
        
        Spectrum contribution = transmittance * geo_term * bsdf * emitted / illumination_pdf;
        Real scatter_pdf = pdf_sample_phase(phase_func, inc_dir, -outgoing_dir) * geo_term * path_prob;
        Real weight = (illumination_pdf * illumination_pdf) / (illumination_pdf * illumination_pdf + scatter_pdf * scatter_pdf);

        return weight * contribution;
    }

    return make_zero_spectrum();
}



Spectrum nee_brdf(Vector3 pos, Vector3 inc_dir, int med_id, int bounce_cnt, Light src_light, const Scene& scn, pcg32_state& rng_state, int src_id, PathVertex v) {

    PointAndNormal light_sample = sample_point_on_light(src_light, pos, 
                                Vector2(next_pcg32_real<Real>(rng_state), next_pcg32_real<Real>(rng_state)),
                                next_pcg32_real<Real>(rng_state), scn);

    Real transmittance = 1;
    int active_medium = med_id;
    int bounces_made = 0;
    Real path_prob = 1;

    Vector3 start_pos = pos;

    while (true) {
        Vector3 ray_dir = normalize(light_sample.position - pos);
        Ray trace_ray = Ray{pos, ray_dir, get_shadow_epsilon(scn),
                            (1 - get_shadow_epsilon(scn)) * distance(light_sample.position, pos)};
        RayDifferential diff_ray = RayDifferential{ Real(0), Real(0) };
        std::optional<PathVertex> hit = intersect(scn, trace_ray, diff_ray);
        Real ray_length = distance(pos, light_sample.position);

        if (hit) {
            ray_length = distance(pos, hit->position);
        }
        if (active_medium != -1) {
            Medium medium_props = scn.media[active_medium];
            Real absorption = get_sigma_a(medium_props, Vector3(1, 2, 3)).x;
            Real scatter = get_sigma_s(medium_props, Vector3(1, 2, 3)).x;
            Real extinction = absorption + scatter;
            Real decay = exp(-extinction * ray_length);
            transmittance *= decay;
            path_prob *= decay;
        }

        if (!hit) break;

        PathVertex vertex = hit.value();
        if (vertex.material_id >= 0) {
            return make_zero_spectrum();
        }

        bounces_made++;
        if (scn.options.max_depth != -1 && bounce_cnt + bounces_made + 1 >= scn.options.max_depth) {
            return make_zero_spectrum();
        }

        active_medium = update_medium_id(vertex, trace_ray, active_medium);
        pos += ray_length * trace_ray.dir;
    }

    if (transmittance > 0) {
        Vector3 outgoing_dir = normalize(start_pos - light_sample.position);
        Real sq_distance = distance_squared(start_pos, light_sample.position);
        Real cos_theta = abs(dot(outgoing_dir, light_sample.normal));
        Real geo_term = cos_theta / sq_distance;

        const Material& mat = scn.materials[v.material_id];

        Spectrum bsdf = eval(mat, inc_dir, -outgoing_dir, v, scn.texture_pool);
        Spectrum emitted = emission(src_light, outgoing_dir, Real(0), light_sample, scn);
        Real illumination_pdf = light_pmf(scn, src_id) * pdf_point_on_light(src_light, light_sample, start_pos, scn);
        
        Spectrum contribution = transmittance * geo_term * bsdf * emitted / illumination_pdf;
        Real scatter_pdf = pdf_sample_bsdf(mat, inc_dir, -outgoing_dir, v, scn.texture_pool) * geo_term * path_prob;
        Real weight = (illumination_pdf * illumination_pdf) / (illumination_pdf * illumination_pdf + scatter_pdf * scatter_pdf);

        return weight * contribution;
    }

    return make_zero_spectrum();
}


// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    Real X = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real Y = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    Vector2 screen_pos(X, Y);

    Ray camera_ray = sample_primary(scene.camera, screen_pos);
    std::optional<PathVertex> isect = intersect(scene, camera_ray);
    if (isect) {
        PathVertex path_vertex = isect.value();
        if (is_light(scene.shapes[path_vertex.shape_id])) {
            Spectrum transmittance = make_const_spectrum(1);
            Spectrum Le = emission(path_vertex, -camera_ray.dir, scene);
            if (scene.camera.medium_id >= 0) {
                Real t = distance(camera_ray.org, path_vertex.position);
                Spectrum sigma_a = get_sigma_a(scene.media[path_vertex.exterior_medium_id], path_vertex.position);
                transmittance = exp(-sigma_a * t);
            }
            return transmittance * Le;
        }
    }
    return make_zero_spectrum();
}


// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    Real X = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real Y = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    Vector2 screen_pos(X, Y);

    Ray camera_ray = sample_primary(scene.camera, screen_pos);
    int current_medium_id = scene.camera.medium_id;
    std::optional<PathVertex> isect = intersect(scene, camera_ray);
    Spectrum radiance = make_zero_spectrum();

    if (current_medium_id >= 0) {
        // PathVertex path_vertex = isect.value();
        const Medium &medium = scene.media[current_medium_id];
        Spectrum sigma_a = get_sigma_a(medium, camera_ray.org);
        Spectrum sigma_s = get_sigma_s(medium, camera_ray.org);
        Spectrum sigma_t = sigma_s + sigma_a;

        Real t_hit = infinity<Real>();
        if (isect) {
            t_hit = distance(isect->position, camera_ray.org);
        }
        Real u = next_pcg32_real<Real>(rng);
        Real t = -log(1 - u) / sigma_t.x;
        Vector3 p = camera_ray.org + t * camera_ray.dir;

        if (t < t_hit) {
            PathVertex vertex;
            vertex.position = p;
            vertex.interior_medium_id = current_medium_id;
            vertex.exterior_medium_id = current_medium_id;

            Spectrum transmittance = Real(1) / sigma_t;
            const PhaseFunction &phase_function = get_phase_function(medium);
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal point_on_light =
                sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);
            Real G = 0;
            Vector3 dir_light = normalize(point_on_light.position - vertex.position);
            Ray shadow_ray{vertex.position, dir_light, 
                            get_shadow_epsilon(scene),
                            (1 - get_shadow_epsilon(scene)) *
                                distance(point_on_light.position, vertex.position)};
            if (!occluded(scene, shadow_ray)) {
                G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                    distance_squared(point_on_light.position, vertex.position);
            }
            Real p1 = light_pmf(scene, light_id) *
                pdf_point_on_light(light, point_on_light, vertex.position, scene);
            if (G > 0 && p1 > 0) {
                Vector3 dir_view = -camera_ray.dir;
                Spectrum f = eval(phase_function, dir_view, dir_light);
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                Spectrum T = exp(-sigma_t * distance(vertex.position, point_on_light.position));
                radiance = sigma_s * transmittance * T * G * f * L / p1;
            }
        } else {
            if (isect) {
                const PathVertex path_vertex = isect.value();
                if (is_light(scene.shapes[path_vertex.shape_id])) {
                    radiance = emission(path_vertex, -camera_ray.dir, scene);
                }
            }
        }
    }
    return radiance;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
  
    // 计算随机采样的屏幕空间坐标
    Real screenX = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real screenY = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    // 从相机生成一条主光线
    Ray primaryRay = sample_primary(scene.camera, Vector2(screenX, screenY));
    // 初始化当前路径的透过率
    Spectrum current_path_throughput(1.0, 1.0, 1.0);
    Spectrum radiance(0.0, 0.0, 0.0);

    int bounces  = 0;
    // 初始化射线的起始介质
    int currentMediumIndex = scene.camera.medium_id;

    while (true) {
        // 是否发生散射
        bool scatter = false; 
        // 检测光线与场景的交点
        std::optional<PathVertex> isect = intersect(scene, primaryRay); // 检测光线与场景的交点,可能不与表面相交，但我们可能在一种介质内（在else里面）
        Real transmittance = 1; //散射事件发生的概率
        Real trans_pdf = 1;  //光线在介质中传播时被吸收的程度
        if (currentMediumIndex != -1) { //如果不是在真空中
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real sigmaS = get_sigma_s(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaA = get_sigma_a(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaT = sigmaA + sigmaS;
            Real randomSample = next_pcg32_real<Real>(rng);
            Real t = -std::log(1.0 - randomSample) / sigmaT;
            // compute transmittance and trans_pdf
            if (isect) { // 如果光线与场景相交
                Real t_hit = distance(isect->position, primaryRay.org); 
                if (t < t_hit) { //光线在到达表面之前就已经散射。
                    scatter = true;
                    trans_pdf = exp(-sigmaT * t) * sigmaT;
                    transmittance = exp(-sigmaT * t);
                    primaryRay.org = primaryRay.org + t * primaryRay.dir;
                } else { //光线首先与表面相交
                    scatter = false;
                    trans_pdf = exp(-sigmaT * t_hit);
                    transmittance = exp(-sigmaT * t_hit);
                    primaryRay.org = primaryRay.org + t_hit * primaryRay.dir;
                } 
            } else { // 光线没有与场景相交,直接在介质中散射
                scatter = true;
                trans_pdf = exp(-sigmaT * t) * sigmaT;
                transmittance = exp(-sigmaT * t);
                primaryRay.org = primaryRay.org + t * primaryRay.dir;
            }  
        }       
            
        if (!scatter) { // 如果没有发生散射，而是到达了一个会发光的表面
            if (isect) { // 如果光线与场景相交
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) {  //如果交点是光源
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene);
                }
            }
        }

        if (scene.options.max_depth > 0 && bounces == scene.options.max_depth-1) break;
        //处理路径的一次bounce
        if (!scatter && isect) { // 如果没有发生散射，而是到达了一个对光线没有产生影响的表面
            PathVertex vertex = isect.value();
            if (vertex.material_id == -1) {
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                primaryRay.org = vertex.position + primaryRay.dir * get_intersection_epsilon(scene);
                bounces += 1;
                continue;
            }
        }

        
        if (scatter) {
            //# sample next direct & update path throughput
            Vector3 p = primaryRay.org;
            std::optional<Vector3> opt_dir = sample_phase_function(get_phase_function(scene.media[currentMediumIndex]), -primaryRay.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
            if (opt_dir) {
                Vector3 next_dir = opt_dir.value();
                Real sigma_s = get_sigma_s(scene.media[currentMediumIndex], p).x;
                Real pdf = pdf_sample_phase(get_phase_function(scene.media[currentMediumIndex]), -primaryRay.dir, next_dir);
                current_path_throughput *= eval(get_phase_function(scene.media[currentMediumIndex]), -primaryRay.dir, next_dir) / pdf;
                primaryRay = Ray{primaryRay.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
            }
        }
        else {
            //Hit a surface -- don’t need to deal with this yet
            break;
        }

        // 检查是否需要进行俄罗斯罗盘
        // 实施俄罗斯轮盘算法
        if (scene.options.rr_depth > 0 && bounces >= scene.options.rr_depth) {
            Real rr_prob = std::min(current_path_throughput.x, 0.95);
            Real randomValue = next_pcg32_real<Real>(rng);
            if (rr_prob < randomValue) break;
            current_path_throughput /= rr_prob;
        }
        bounces++;
    }
    return radiance;
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    
    // 计算随机采样的屏幕空间坐标
    Real screenX = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real screenY = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    // 从相机生成一条主光线
    Ray primaryRay = sample_primary(scene.camera, Vector2(screenX, screenY));
    // 初始化当前路径的透过率
    Spectrum current_path_throughput(1.0, 1.0, 1.0);
    Spectrum radiance(0.0, 0.0, 0.0);

    int bounces  = 0;
    // 初始化射线的起始介质
    int currentMediumIndex = scene.camera.medium_id;

    bool never_scatter = true;
    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1;

    while (true) {
        // 是否发生散射
        bool scatter = false; 
        // 检测光线与场景的交点
        std::optional<PathVertex> isect = intersect(scene, primaryRay); // 检测光线与场景的交点,可能不与表面相交，但我们可能在一种介质内（在else里面）
        Real transmittance = 1; //散射事件发生的概率
        Real trans_pdf = 1;  //光线在介质中传播时被吸收的程度
        
        if (currentMediumIndex == -1) { // 在真空中，直接处理
            if (isect) {
                PathVertex vertex = isect.value();
                primaryRay.org = primaryRay.org + distance(primaryRay.org, vertex.position) * primaryRay.dir;
            } else break;
        }
        
        if (currentMediumIndex != -1) { // 不在真空中，处理在介质中的散射、吸收情况
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real sigmaS = get_sigma_s(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaA = get_sigma_a(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaT = sigmaA + sigmaS;
            Real randomSample = next_pcg32_real<Real>(rng);
            Real t = -std::log(1.0 - randomSample) / sigmaT; // 使用指数分布的逆变换采样 t
            // compute transmittance and trans_pdf
            if (isect) { // 如果光线与场景相交
                Real t_hit = distance(isect->position, primaryRay.org); 
                if (t < t_hit) { //光线在到达表面之前就已经散射。
                    scatter = true;
                    trans_pdf = exp(-sigmaT * t) * sigmaT;
                    transmittance = exp(-sigmaT * t);
                    primaryRay.org = primaryRay.org + t * primaryRay.dir;
                } else { //光线首先与表面相交
                    scatter = false;
                    trans_pdf = exp(-sigmaT * t_hit);
                    transmittance = exp(-sigmaT * t_hit);
                    primaryRay.org = primaryRay.org + t_hit * primaryRay.dir;
                } 
            } else { // 光线没有与场景相交,直接在介质中散射
                scatter = true;
                trans_pdf = exp(-sigmaT * t) * sigmaT;
                transmittance = exp(-sigmaT * t);
                primaryRay.org = primaryRay.org + t * primaryRay.dir;
            }  
        }       
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && isect) { // 没有散射而到达一个光源
            if (never_scatter) {
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) {  // 如果交点是光源
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene);
                }
            } else { // 之前有至少发生一次过散射
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) { // 多重采样
                    PointAndNormal light_point = {vertex.position, vertex.geometric_normal};
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light currlight = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(currlight, light_point, nee_p_cache, scene) * light_pmf(scene, light_id); // 采样到此光源的总概率
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache); // 最后一次散射位置 nee_p_cache 到光源的向量
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(vertex.position - nee_p_cache);
                    Real G = top/bottom;
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene) * w;
                }

            }
        }

        if (scene.options.max_depth > 0 && bounces == scene.options.max_depth-1) break;
        //处理路径的一次bounce
        if (!scatter && isect) {
            PathVertex vertex = isect.value();
            if (vertex.material_id == -1) { // 如果没有发生散射，而是到达了一个对光线没有产生影响的表面
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                primaryRay.org = vertex.position + primaryRay.dir * get_intersection_epsilon(scene);
                bounces += 1;
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        
        if (scatter) {
            never_scatter = false;
            // sample next direct & update path throughput
            Vector3 p = primaryRay.org;
            PhaseFunction pf = get_phase_function(scene.media[currentMediumIndex]);
            std::optional<Vector3> opt_dir = sample_phase_function(pf, -primaryRay.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
            if (opt_dir) {
                Vector3 next_dir = opt_dir.value();
                Real sigma_s = get_sigma_s(scene.media[currentMediumIndex], p).x;
                int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                Spectrum nee_out = next_event_estimation(p, -primaryRay.dir, currentMediumIndex, bounces, scene.lights[light_id], scene, rng, light_id);
                radiance += current_path_throughput * nee_out * sigma_s;
                dir_pdf = pdf_sample_phase(pf, -primaryRay.dir, next_dir);
                current_path_throughput *= (eval(pf, -primaryRay.dir, next_dir) / dir_pdf) * sigma_s;
                primaryRay = Ray{primaryRay.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                nee_p_cache = p;;
                multi_trans_pdf = Real(1);
            }
        }
        else {
            //Hit a surface -- don’t need to deal with this yet
            break;
        }

        // 检查是否需要进行俄罗斯罗盘
        // 实施俄罗斯轮盘算法
        if (scene.options.rr_depth > 0 && bounces >= scene.options.rr_depth) {
            Real rr_prob = std::min(current_path_throughput.x, 0.95);
            Real randomValue = next_pcg32_real<Real>(rng);
            if (rr_prob < randomValue) break;
            current_path_throughput /= rr_prob;
        }
        bounces++;
    }
    return radiance;
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // 计算随机采样的屏幕空间坐标
    Real screenX = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real screenY = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    // 从相机生成一条主光线
    Ray primaryRay = sample_primary(scene.camera, Vector2(screenX, screenY));
    // 初始化当前路径的透过率
    Spectrum current_path_throughput(1.0, 1.0, 1.0);
    Spectrum radiance(0.0, 0.0, 0.0);

    int bounces  = 0;
    // 初始化射线的起始介质
    int currentMediumIndex = scene.camera.medium_id;

    bool never_scatter = true;
    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Real multi_trans_pdf = 1;

    while (true) {
        // 是否发生散射
        bool scatter = false; 
        // 检测光线与场景的交点
        std::optional<PathVertex> isect = intersect(scene, primaryRay); // 检测光线与场景的交点,可能不与表面相交，但我们可能在一种介质内（在else里面）
        Real transmittance = 1; //散射事件发生的概率
        Real trans_pdf = 1;  //光线在介质中传播时被吸收的程度
        
        if (currentMediumIndex == -1) { // 在真空中，直接处理
            if (isect) {
                PathVertex vertex = isect.value();
                primaryRay.org = primaryRay.org + distance(primaryRay.org, vertex.position) * primaryRay.dir;
            } else break;
        }
        
        if (currentMediumIndex != -1) { // 不在真空中，处理在介质中的散射、吸收情况
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Real sigmaS = get_sigma_s(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaA = get_sigma_a(scene.media[currentMediumIndex], primaryRay.org).x;
            Real sigmaT = sigmaA + sigmaS;
            Real randomSample = next_pcg32_real<Real>(rng);
            Real t = -std::log(1.0 - randomSample) / sigmaT; // 使用指数分布的逆变换采样 t
            // compute transmittance and trans_pdf
            if (isect) { // 如果光线与场景相交
                Real t_hit = distance(isect->position, primaryRay.org); 
                if (t < t_hit) { //光线在到达表面之前就已经散射。
                    scatter = true;
                    trans_pdf = exp(-sigmaT * t) * sigmaT;
                    transmittance = exp(-sigmaT * t);
                    primaryRay.org = primaryRay.org + t * primaryRay.dir;
                } else { //光线首先与表面相交
                    scatter = false;
                    trans_pdf = exp(-sigmaT * t_hit);
                    transmittance = exp(-sigmaT * t_hit);
                    primaryRay.org = primaryRay.org + t_hit * primaryRay.dir;
                } 
            } else { // 光线没有与场景相交,直接在介质中散射
                scatter = true;
                trans_pdf = exp(-sigmaT * t) * sigmaT;
                transmittance = exp(-sigmaT * t);
                primaryRay.org = primaryRay.org + t * primaryRay.dir;
            }  
        }       
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && isect) { // 没有散射而到达一个光源
            if (never_scatter) {
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) {  // 如果交点是光源
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene);
                }
            } else { // 之前有至少发生一次过散射
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) { // 多重采样
                    PointAndNormal light_point = {vertex.position, vertex.geometric_normal};
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light currlight = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(currlight, light_point, nee_p_cache, scene) * light_pmf(scene, light_id); // 采样到此光源的总概率
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache); // 最后一次散射位置 nee_p_cache 到光源的向量
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(vertex.position - nee_p_cache);
                    Real G = top/bottom;
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene) * w;
                }

            }
        }

        if (scene.options.max_depth > 0 && bounces == scene.options.max_depth-1) break;
        //处理路径的一次bounce
        if (!scatter && isect) {
            PathVertex vertex = isect.value();
            if (vertex.material_id == -1) { // 如果没有发生散射，而是到达了一个对光线没有产生影响的表面
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                primaryRay.org = vertex.position + primaryRay.dir * get_intersection_epsilon(scene);
                bounces += 1;
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        
        if (scatter) {
            never_scatter = false;
            // sample next direct & update path throughput
            Vector3 p = primaryRay.org;
            PhaseFunction pf = get_phase_function(scene.media[currentMediumIndex]);
            std::optional<Vector3> opt_dir = sample_phase_function(pf, -primaryRay.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
            if (opt_dir) {
                Vector3 next_dir = opt_dir.value();
                Real sigma_s = get_sigma_s(scene.media[currentMediumIndex], p).x;
                int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                Spectrum nee_out = next_event_estimation(p, -primaryRay.dir, currentMediumIndex, bounces, scene.lights[light_id], scene, rng, light_id);
                radiance += current_path_throughput * nee_out * sigma_s;
                dir_pdf = pdf_sample_phase(pf, -primaryRay.dir, next_dir);
                current_path_throughput *= (eval(pf, -primaryRay.dir, next_dir) / dir_pdf) * sigma_s;
                primaryRay = Ray{primaryRay.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                nee_p_cache = p;;
                multi_trans_pdf = Real(1);
            }
        }
        else {

            if (isect) {
                PathVertex vertex = isect.value();

                int sampled_light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                Spectrum light_contribution = nee_brdf(primaryRay.org, -primaryRay.dir, currentMediumIndex, bounces, scene.lights[sampled_light_id], scene, rng, sampled_light_id, vertex);
                radiance += current_path_throughput * light_contribution;
                never_scatter = false;

                const Material &material = scene.materials[vertex.material_id];
                Vector3 view_dir = -primaryRay.dir;
                Vector2 bsdf_sample_params{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_sample_weight = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_result = sample_bsdf(material, view_dir, vertex, scene.texture_pool, bsdf_sample_params, bsdf_sample_weight);

                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_result;
                Vector3 outgoing_dir = bsdf_sample.dir_out;
                primaryRay = Ray{primaryRay.org, outgoing_dir, get_intersection_epsilon(scene), infinity<Real>()};

                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);

                Spectrum bsdf_value = eval(material, view_dir, outgoing_dir, vertex, scene.texture_pool);
                Real bsdf_pdf_value = pdf_sample_bsdf(material, view_dir, outgoing_dir, vertex, scene.texture_pool);
                current_path_throughput *= bsdf_value / bsdf_pdf_value;

                // Update cache with new values
                dir_pdf = bsdf_pdf_value;
                nee_p_cache = primaryRay.org;
                multi_trans_pdf = Real(1);
            }

        }

        // 检查是否需要进行俄罗斯罗盘
        // 实施俄罗斯轮盘算法
        if (scene.options.rr_depth > 0 && bounces >= scene.options.rr_depth) {
            Real rr_prob = std::min(current_path_throughput.x, 0.95);
            Real randomValue = next_pcg32_real<Real>(rng);
            if (rr_prob < randomValue) break;
            current_path_throughput /= rr_prob;
        }
        bounces++;
    }
    return radiance;
}


const Real r_weight = 0.11060479323; //for 680nm 
const Real g_weight = 0.25843946974; //for 550nm
const Real b_weight = 0.63095573702; //for 440nm

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    // Homework 2: implememt this!

    // 计算随机采样的屏幕空间坐标
    Real screenX = (next_pcg32_real<Real>(rng) + x) / scene.camera.width;
    Real screenY = (next_pcg32_real<Real>(rng) + y) / scene.camera.height;
    // 从相机生成一条主光线
    Ray primaryRay = sample_primary(scene.camera, Vector2(screenX, screenY));
    // 初始化当前路径的透过率
    Spectrum current_path_throughput(1.0, 1.0, 1.0);
    Spectrum radiance(0.0, 0.0, 0.0);
    
    int bounces  = 0;
    // 初始化射线的起始介质
    int currentMediumIndex = scene.camera.medium_id;

    bool never_scatter = true;
    Real dir_pdf = 0;
    Vector3 nee_p_cache;
    Spectrum multi_trans_pdf = make_const_spectrum(1);

    int max_depth = scene.options.max_depth;

    for (bounces = 0; max_depth == -1 || bounces < max_depth; bounces++) {
        // 是否发生散射
        bool scatter = false; 
        // 检测光线与场景的交点
        std::optional<PathVertex> isect = intersect(scene, primaryRay); // 检测光线与场景的交点,可能不与表面相交，但我们可能在一种介质内（在else里面）
        Spectrum transmittance = make_const_spectrum(1); //散射事件发生的概率
        Spectrum trans_pdf = make_const_spectrum(1); //光线在介质中传播时被吸收的程度
        
        PathVertex vertex;
        if (isect) {
            vertex = *isect;
        }

        if (currentMediumIndex == -1) { // 在真空中，直接处理
            if (isect) {
                PathVertex vertex = isect.value();
                primaryRay.org = primaryRay.org + distance(primaryRay.org, vertex.position) * primaryRay.dir;
            } else break;
        }
        
        if (currentMediumIndex != -1) { // 不在真空中，处理在介质中的散射、吸收情况
            // sample t s.t. p(t) ~ exp(-sigma_t * t)
            Spectrum sigma_a = get_sigma_a(scene.media[currentMediumIndex], primaryRay.org);
            Spectrum sigma_s = get_sigma_s(scene.media[currentMediumIndex], primaryRay.org);
            Spectrum sigma_t = sigma_s + sigma_a;
      
            // Importance Sample a channel for sampling
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * Real(3)), 0, 2);
            // do not Assume monochromatic medium
            Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[channel];

            if (isect) { // 如果光线与场景相交
                Real t_hit = distance(isect->position, primaryRay.org); 
                if (t < t_hit) { //光线在到达表面之前就已经散射。
                    scatter = true;
                    never_scatter = false;
                    trans_pdf = exp(-sigma_t * t) * sigma_t;
                    transmittance = exp(-sigma_t * t);

                    // Update vertex information
                    vertex.position = primaryRay.org + t * primaryRay.dir;
                    vertex.interior_medium_id = currentMediumIndex;
                    vertex.exterior_medium_id = currentMediumIndex;
                    // primaryRay.org = primaryRay.org + t * primaryRay.dir;
                } else { //光线首先与表面相交
                    scatter = false;
                    trans_pdf = exp(-sigma_t * t_hit);
                    transmittance = exp(-sigma_t * t_hit);
                    // primaryRay.org = primaryRay.org + t_hit * primaryRay.dir;
                } 
            } else { // 光线没有与场景相交,直接在介质中散射
                // scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                // primaryRay.org = primaryRay.org + t * primaryRay.dir;
            }  
        }
        multi_trans_pdf *= trans_pdf;
        current_path_throughput *= (transmittance / average(trans_pdf));

        if (!scatter && !isect) {
            break;
        }

        if (!scatter && isect) { // 没有散射而到达一个光源
            if (never_scatter) {
                // PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) { // 如果交点是光源
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene);
                }
            } else {
                // PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) { // 多重采样
                    PointAndNormal light_point{vertex.position, vertex.geometric_normal};
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light currlight = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(currlight, light_point, nee_p_cache, scene) * light_pmf(scene, light_id); // 采样到此光源的总概率
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache); // 最后一次散射位置 nee_p_cache 到光源的向量
                    Real top = abs(dot(omega_prime, vertex.geometric_normal));
                    Real bottom = length_squared(vertex.position - nee_p_cache);
                    Real G = top/bottom;
                    Real p_dir = dir_pdf * average(multi_trans_pdf) * G;
                    Real w = (p_dir * p_dir) / (p_dir * p_dir + pdf_nee * pdf_nee);
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene) * w;
                }
                
            }
        }

        if (scene.options.max_depth > 0 && bounces == scene.options.max_depth-1) break;
        //处理路径的一次bounce
        if (!scatter && isect) {
            PathVertex vertex = isect.value();
            if (vertex.material_id == -1) { // 如果没有发生散射，而是到达了一个对光线没有产生影响的表面
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                primaryRay.org = vertex.position + primaryRay.dir * get_intersection_epsilon(scene);
                // bounces += 1;
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }


        // Tracer 6 and atmosphere

        // Cache the NEE vertex for later use
        nee_p_cache = vertex.position;
        // We have a scattering event (medium or surface), reset multi_trans_pdf
        multi_trans_pdf = make_const_spectrum(1);

        // next event estimation
        Spectrum C1 = make_zero_spectrum();
        Real w1 = 0;
        {
            // Sample a point on the light source.
            Vector2 light_uv{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
            Real light_w = next_pcg32_real<Real>(rng);
            Real shape_w = next_pcg32_real<Real>(rng);
            int light_id = sample_light(scene, light_w);
            const Light &light = scene.lights[light_id];
            PointAndNormal point_on_light =
                sample_point_on_light(light, vertex.position, light_uv, shape_w, scene);
            // Compute transmittance to light. Skip through index-matching shapes.
            Spectrum T_light = make_const_spectrum(1);
            Vector3 p = vertex.position;
            int shadow_medium_id = currentMediumIndex;
            int shadow_bounces = 0;
            Spectrum p_trans_dir = make_const_spectrum(1);
            while (true) {
                Vector3 dir_light = normalize(point_on_light.position - p);
                Ray shadow_ray{p, dir_light, 
                               get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(point_on_light.position, p)};
                std::optional<PathVertex> shadow_vertex = intersect(scene, shadow_ray);
                Real next_t = shadow_ray.tfar;
                if (shadow_vertex) {
                    next_t = distance(p, shadow_vertex->position);
                }

                // Account for the transmittance to next_t
                if (shadow_medium_id >= 0) {
                    const Medium &medium = scene.media[shadow_medium_id];
                    Spectrum sigma_a = get_sigma_a(medium, primaryRay.org);
                    Spectrum sigma_s = get_sigma_s(medium, primaryRay.org);
                    Spectrum sigma_t = sigma_s + sigma_a;

                    T_light *= exp(-sigma_t * next_t);
                    p_trans_dir *= exp(-sigma_t * next_t);
                }

                if (!shadow_vertex) {
                    // Nothing is blocking, we're done
                    break;
                } else {
                    // Something is blocking: is it an opaque surface?
                    if (shadow_vertex->material_id >= 0) {
                        // we're blocked
                        T_light = make_zero_spectrum();
                        p_trans_dir = make_zero_spectrum();
                        break;
                    }
                    // otherwise, we want to pass through -- this introduces
                    // one extra connection vertex
                    shadow_bounces++;
                    if (max_depth != -1 && bounces + shadow_bounces + 1 >= max_depth) {
                        // Reach the max no. of vertices
                        T_light = make_zero_spectrum();
                        break;
                    }
                    // let's update and continue
                    shadow_medium_id = update_medium_id(*shadow_vertex, shadow_ray, shadow_medium_id);
                    p = p + next_t * dir_light;
                }
            }
            
            if (max(T_light) > 0) {
                // Compute sigma_s * T * T_light * G * f * L
                Vector3 dir_light = normalize(point_on_light.position - vertex.position);
                Real G = max(-dot(dir_light, point_on_light.normal), Real(0)) /
                        distance_squared(point_on_light.position, vertex.position);
                Real p1 = light_pmf(scene, light_id) *
                    pdf_point_on_light(light, point_on_light, vertex.position, scene);
                Vector3 dir_view = -primaryRay.dir;
                Spectrum f;
                // are we on a surface or are we in a medium?
                if (scatter) {
                    assert(currentMediumIndex >= 0);
                    const Medium &medium = scene.media[currentMediumIndex];
                    const PhaseFunction &phase_function = get_phase_function(medium);
                    f = eval(phase_function, dir_view, dir_light);
                } else {
                    const Material &mat = scene.materials[vertex.material_id];
                    f = eval(mat, dir_view, dir_light, vertex, scene.texture_pool);
                }
                Spectrum sigma_s = make_const_spectrum(1);
                if (scatter) {
                    const Medium &medium = scene.media[currentMediumIndex];
                    sigma_s = get_sigma_s(medium, vertex.position);
                }
                Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                C1 = current_path_throughput * sigma_s * T_light * G * f * L / p1;
                // Multiple importance sampling: it's also possible
                // that a phase function sampling + multiple steps exponential sampling
                // will reach the light source.
                // The probability for multiple steps exponential sampling
                // is stored in p_trans_dir
                // We also need to multiply with G to convert phase function PDF to area measure.
                Real p2 = 0;
                if (scatter) {
                    assert(currentMediumIndex >= 0);
                    const Medium &medium = scene.media[currentMediumIndex];
                    const PhaseFunction &phase_function = get_phase_function(medium);
                    p2 = pdf_sample_phase(phase_function, dir_view, dir_light) * G;
                } else {
                    assert(vertex.material_id >= 0);
                    const Material &mat = scene.materials[vertex.material_id];
                    p2 = pdf_sample_bsdf(mat, dir_view, dir_light, vertex, scene.texture_pool) * G;
                }
                p2 *= average(p_trans_dir);
                w1 = (p1 * p1) / (p1 * p1 + p2 * p2);
            }
        }
        radiance += C1 * w1;



        Vector3 next_dir;
        if (scatter) {
            // never_scatter = false;
            // sample next direct & update path throughput
            Vector3 p = primaryRay.org;
            PhaseFunction pf = get_phase_function(scene.media[currentMediumIndex]);
            std::optional<Vector3> opt_dir = sample_phase_function(pf, -primaryRay.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng));
            if (opt_dir) {
                next_dir = opt_dir.value();
                Spectrum sigma_s = get_sigma_s(scene.media[currentMediumIndex], vertex.position);
                // int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                // Spectrum nee_out = next_event_estimation(p, -primaryRay.dir, currentMediumIndex, bounces, scene.lights[light_id], scene, rng, light_id);
                // radiance += current_path_throughput * nee_out * sigma_s;
                dir_pdf = pdf_sample_phase(pf, -primaryRay.dir, next_dir);
                current_path_throughput *= (eval(pf, -primaryRay.dir, next_dir) / dir_pdf) * sigma_s;
                // primaryRay = Ray{primaryRay.org + next_dir * get_intersection_epsilon(scene), next_dir, Real(0), infinity<Real>()};
                // nee_p_cache = p;;
                // multi_trans_pdf = Real(1);
            }
        } 
        else {
            if (isect) {
                // PathVertex vertex = isect.value();

                // int sampled_light_id = sample_light(scene, next_pcg32_real<Real>(rng));
                // Spectrum light_contribution = nee_brdf(primaryRay.org, -primaryRay.dir, currentMediumIndex, bounces, scene.lights[sampled_light_id], scene, rng, sampled_light_id, vertex);
                // radiance += current_path_throughput * light_contribution;
                never_scatter = false;

                const Material &material = scene.materials[vertex.material_id];
                Vector3 view_dir = -primaryRay.dir;
                Vector2 bsdf_sample_params{next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)};
                Real bsdf_sample_weight = next_pcg32_real<Real>(rng);
                std::optional<BSDFSampleRecord> bsdf_sample_result = sample_bsdf(material, view_dir, vertex, scene.texture_pool, bsdf_sample_params, bsdf_sample_weight);

                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_result;
                Vector3 outgoing_dir = bsdf_sample.dir_out;
                // primaryRay = Ray{primaryRay.org, outgoing_dir, get_intersection_epsilon(scene), infinity<Real>()};
                
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                
                Spectrum bsdf_value = eval(material, view_dir, outgoing_dir, vertex, scene.texture_pool);
                Real bsdf_pdf_value = pdf_sample_bsdf(material, view_dir, outgoing_dir, vertex, scene.texture_pool);
                current_path_throughput *= bsdf_value / bsdf_pdf_value;
                
                // Update cache with new values
                // dir_pdf = bsdf_pdf_value;
                // nee_p_cache = primaryRay.org;
                // multi_trans_pdf = Real(1);
            }
        }

        // 检查是否需要进行俄罗斯罗盘
        // 实施俄罗斯轮盘算法
        if (scene.options.rr_depth > 0 && bounces >= scene.options.rr_depth) {
            Real rr_prob = std::min(current_path_throughput.x, 0.95);
            Real randomValue = next_pcg32_real<Real>(rng);
            if (rr_prob < randomValue) break;
            current_path_throughput /= rr_prob;
        }

        // Update rays
        primaryRay = Ray{vertex.position,
                  next_dir,
                  get_intersection_epsilon(scene),
                  infinity<Real>()};
        currentMediumIndex =
            update_medium_id(vertex, primaryRay, currentMediumIndex);
    }
    return radiance;
}
