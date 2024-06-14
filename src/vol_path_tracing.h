#pragma once

#include "lajolla.h"
#include "pcg.h"
#include "phase_function.h"
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

Spectrum next_event_estimation(Vector3 p, Vector3 omega, int current_medium, int bounces, Light light, const Scene& scene, pcg32_state& rng, int light_id) {

    PointAndNormal p_prime = sample_point_on_light(light, p, 
                                Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                next_pcg32_real<Real>(rng), scene);

    Real T_light = 1;
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1;

    Vector3 orig_p = p;

    while (true) {
        Ray shadow_ray = Ray{p, normalize(p_prime.position - p),  get_shadow_epsilon(scene),
                               (1 - get_shadow_epsilon(scene)) *
                                   distance(p_prime.position, p)};
        RayDifferential ray_diff = RayDifferential{ Real(0), Real(0) };
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);

        if (isect) {
            PathVertex vertex = *isect;
            next_t = distance(p, vertex.position);
        }

        if (shadow_medium != -1) {
            Medium medium = scene.media[shadow_medium];
            Real sigma_a = get_sigma_a(medium, Vector3(1, 2, 3)).x;
            Real sigma_s = get_sigma_s(medium, Vector3(1, 2, 3)).x;
            Real sigma_t = sigma_a + sigma_s;
            T_light *= exp(-sigma_t * next_t);
            p_trans_dir *= exp(-sigma_t * next_t);
        }
            
        if (isect) {
            PathVertex vertex = isect.value();
            if (vertex.material_id >= 0)
                return make_zero_spectrum();
            ++shadow_bounces;
            int max_bounce = scene.options.max_depth;
            if (max_bounce != -1 && bounces + shadow_bounces + 1 >= max_bounce)
                return make_zero_spectrum();
            shadow_medium = update_medium_id(vertex, shadow_ray, 0);
            p = p + next_t * (shadow_ray.dir);
        } else break;
    }

    if (T_light > 0) {
        Vector3 omega_prime = normalize(orig_p - p_prime.position);
        Real denom = distance_squared(orig_p, p_prime.position);
        Real top = abs(dot(omega_prime, p_prime.normal));
        Real G = top / denom;
        PhaseFunction pf = get_phase_function(scene.media[current_medium]);

        Spectrum f = eval(pf, omega, -omega_prime);
        Spectrum Le = emission(light, omega_prime, Real(0), p_prime, scene);
        Real pdf_nee = light_pmf(scene, light_id) *
            pdf_point_on_light(light, p_prime, orig_p, scene);
        
        Spectrum contrib = T_light * G * f * Le / pdf_nee;
        Real pdf_phase = pdf_sample_phase(pf, omega, -omega_prime) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);

        return w * contrib;
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
        
        if (currentMediumIndex != -1) { // 不在真空中
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
        current_path_throughput *= (transmittance / trans_pdf);

        if (!scatter && isect) { // 如果没有发生散射，而是到达了一个会发光的表面
            if (never_scatter) {
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) {  //如果交点是光源
                    radiance += current_path_throughput * emission(vertex, -primaryRay.dir, scene);
                }
            } else { // 之前有至少发生一次过散射
                PathVertex vertex = isect.value();
                if (is_light(scene.shapes[vertex.shape_id])) {
                    PointAndNormal light_point = {vertex.position, vertex.geometric_normal};
                    int light_id = get_area_light_id(scene.shapes[vertex.shape_id]);
                    Light currlight = scene.lights[light_id];
                    Real pdf_nee = pdf_point_on_light(currlight, light_point, nee_p_cache, scene) * light_pmf(scene, light_id);
                    Vector3 omega_prime = normalize(vertex.position - nee_p_cache);
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
        if (!scatter && isect) { // 如果没有发生散射，而是到达了一个对光线没有产生影响的表面
            PathVertex vertex = isect.value();
            if (vertex.material_id == -1) {
                currentMediumIndex = update_medium_id(vertex, primaryRay, currentMediumIndex);
                primaryRay.org = vertex.position + primaryRay.dir * get_intersection_epsilon(scene);
                bounces += 1;
                multi_trans_pdf *= trans_pdf;
                continue;
            }
        }

        
        if (scatter) {
            never_scatter = false;
            //# sample next direct & update path throughput
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
    return make_zero_spectrum();
}

// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    return make_zero_spectrum();
}
