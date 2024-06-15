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
