## ExoplanetsSysSim/src/eval_model_macro.jl
## (c) 2015 Eric B. Ford

# Macro to create functions that access variables with global scope (i.e., within the current module), so they can be used inside the evaluate model functions
macro make_evaluate_model(param,cat_phys,cat_obs,sum_stat,sum_stat_ref)
  @eval begin
    function evaluate_model(param_vector::Vector{Float64})
      global $param
      update_sim_param_from_vector!(param_vector,$param)
      global $cat_phys = generate_kepler_physical_catalog($param)
      #global $cat_obs = observe_kepler_targets_sky_avg($cat_phys,$param)
      global $cat_obs = observe_kepler_targets_single_obs($cat_phys,$param)
      global $sum_stat = calc_summary_stats_sim_pass_one_demo($cat_obs,$cat_phys,$param)
      dist1 = calc_distance_vector_demo($sum_stat_ref,$sum_stat, 1, $param)
      #if haskey($param,"minimum ABC dist skip pass 2") 
      #   if calc_scalar_distance(dist1) > get($param,"minimum ABC dist skip pass 2") 
      #      return dist1
      #   end
      #end
      $sum_stat = calc_summary_stats_sim_pass_two_demo($cat_obs,$cat_phys,$sum_stat,$param)
      dist2 = calc_distance_vector_demo($sum_stat_ref,$sum_stat, 2, $param)
      return [dist1; dist2]
    end

    function evaluate_model_pass_one(param_vector::Vector{Float64})
      global $param
      update_sim_param_from_vector!(param_vector,$param)
      global $cat_phys = generate_kepler_physical_catalog($param)
      #global $cat_obs = observe_kepler_targets_sky_avg($cat_phys,$param)
      global $cat_obs = observe_kepler_targets_single_obs($cat_phys,$param)
      global $sum_stat = calc_summary_stats_sim_pass_one_demo($cat_obs,$cat_phys,$param)
      dist1 = calc_distance_vector_demo($sum_stat_ref,$sum_stat, 1, $param)
    end

    function evaluate_model_pass_two(param_vector::Vector{Float64})
      global $param
      global $cat_phys
      global $cat_obs
      global $sum_stat = calc_summary_stats_sim_pass_two_demo($cat_obs,$cat_phys,$sum_stat,$param)
      dist2 = calc_distance_vector_demo($sum_stat_ref,$sum_stat, 2, $param)
    end

    function evaluate_model_scalar_ret(param::Vector{Float64})
      calc_scalar_distance(evaluate_model(param))
    end

    function evaluate_model_pass_one_scalar_ret(param::Vector{Float64})
      calc_scalar_distance(evaluate_model_pass_one(param))
    end

    function evaluate_model_pass_two_scalar_ret(param::Vector{Float64})
      calc_scalar_distance(evaluate_model_pass_two(param))
    end

    #if module_name(current_module()) != :Main
    if !isequal(@__MODULE__,Module(:Main))
      export evaluate_model, evaluate_model_pass_one, evaluate_model_pass_two, evaluate_model_scalar_ret, evaluate_model_pass_one_scalar_ret, evaluate_model_pass_two_scalar_ret
    end

  end
end
