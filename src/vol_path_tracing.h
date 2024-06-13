#pragma once
inline int update_medium_id(const PathVertex &vertex,
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
    // Homework 2: implememt this!
    return make_zero_spectrum();
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
