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
int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium_id = scene.camera.medium_id;
    Spectrum radiance = make_zero_spectrum();
    std::optional<PathVertex> vertex_ = intersect(scene, ray, ray_diff);
    if (current_medium_id >= 0) {
        const Medium &medium = scene.media[current_medium_id];
        // We have an integral \int_{0}^{t} T(t') sigma_s f G L dt' + T(max_t) Le
        // T(t) = exp(-\int_{0}^{t} sigma_t dt'')
        // We'll importance sample T
        Real max_t = infinity<Real>();
        if (vertex_) {
            max_t = distance(vertex_->position, ray.org);
        }
        Spectrum sigma_a = get_sigma_a(medium, ray.org);
        Spectrum sigma_s = get_sigma_s(medium, ray.org);
        Spectrum sigma_t = sigma_s + sigma_a;

        // Sample T
        // We want to sample t s.t. p(t) ~ exp(-s * t)
        // We'll do the standard inverse transformation
        // first we integrate from 0 to t
        //   \int_{0}^{t} exp(-s * t') dt'
        // = -exp(-s * t) / s + 1/s
        // the normalization factor is thus 1/s (set t = infty)
        // p(t) = exp(-s * t) * s
        // the CDF is
        // P(t) = -exp(-s * t) + 1 = u
        // to invert the CDF:
        // exp(-s * t) = 1 - u
        // -m * t = log(1 - u)
        // t = log(1 - u) / -s

        // Assume monochromatic medium
        Real t = -log(1 - next_pcg32_real<Real>(rng)) / sigma_t[0];
        if (t < max_t) {
            // Direct lighting
            PathVertex vertex;
            vertex.position = ray.org + t * ray.dir;
            vertex.interior_medium_id = current_medium_id;
            vertex.exterior_medium_id = current_medium_id;

            // Spectrum transmittance = exp(-sigma_t * t);
            // Real pdf = exp(-sigma_t[0] * t) * sigma_t[0];
            // transmittance /= pdf;
            Spectrum transmittance = Real(1) / sigma_t;
            const PhaseFunction &phase_function = get_phase_function(medium);

            // next event estimation
            Spectrum C1 = make_zero_spectrum();
            {
                // Sample a point on the light source.
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
                    // Compute sigma_s * T * G * f * L
                    Vector3 dir_view = -ray.dir;
                    Spectrum f = eval(phase_function, dir_view, dir_light);
                    Spectrum L = emission(light, -dir_light, Real(0), point_on_light, scene);
                    Spectrum T = exp(-sigma_t * distance(vertex.position, point_on_light.position));
                    C1 = sigma_s * transmittance * T * G * f * L / p1;
                }
            }
            radiance += C1;
        } else {
            // Spectrum transmittance = exp(-sigma_t * max_t);
            // Real pdf = exp(-sigma_t[0] * max_t);
            // transmittance /= pdf;
            if (vertex_) {
                const PathVertex vertex = *vertex_;
                if (is_light(scene.shapes[vertex.shape_id])) {
                    radiance += emission(vertex, -ray.dir, scene);
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
