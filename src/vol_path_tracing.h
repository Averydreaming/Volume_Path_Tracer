#pragma once

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
    

    return make_zero_spectrum();
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
