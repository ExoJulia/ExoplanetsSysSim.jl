## ExoplanetsSysSim/src/eval_model.jl
## (c) 2015 Eric B. Ford

# This merely demonstrates how to setup your own function to evaluate a simulated model against an observed catalog
# This code need not be used for any actual calculations
module TestEvalModel
export test_eval_model, evaluate_model_scalar_ret
using ExoplanetsSysSim

sim_param_closure = SimParam()
cat_phys_try_closure = KeplerPhysicalCatalog([])
cat_obs_try_closure = KeplerObsCatalog([])
summary_stat_try_closure =  CatalogSummaryStatistics()
summary_stat_ref_closure =  CatalogSummaryStatistics()

include("eval_model_macro.jl")
@make_evaluate_model(sim_param_closure,cat_phys_try_closure,cat_obs_try_closure,summary_stat_try_closure,summary_stat_ref_closure)

function test_eval_model()
  global sim_param_closure = setup_sim_param_demo()
  cat_phys = generate_kepler_physical_catalog(sim_param_closure)
  cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
  global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
  global cat_phys_try_closure  = generate_kepler_physical_catalog(sim_param_closure)
  global cat_obs_try_closure  = observe_kepler_targets_single_obs(cat_phys_try_closure,sim_param_closure)
  global summary_stat_try_closure  = calc_summary_stats_sim_pass_one_demo(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
  summary_stat_try_closure   = calc_summary_stats_sim_pass_two_demo(cat_obs_try_closure,cat_phys_try_closure,summary_stat_try_closure,sim_param_closure)
  param_guess = make_vector_of_sim_param(sim_param_closure)
  evaluate_model_scalar_ret( param_guess)
end

function test_eval_model_vs_sim_data()
  global sim_param_closure = setup_sim_param_demo()
  cat_phys = generate_kepler_physical_catalog(sim_param_closure)
  # cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
  cat_obs = simulated_read_kepler_observations(sim_param_closure)
  global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
  global cat_phys_try_closure  = generate_kepler_physical_catalog(sim_param_closure)
  global cat_obs_try_closure  = observe_kepler_targets_single_obs(cat_phys_try_closure,sim_param_closure)
  global summary_stat_try_closure  = calc_summary_stats_sim_pass_one_demo(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
  summary_stat_try_closure   = calc_summary_stats_sim_pass_two_demo(cat_obs_try_closure,cat_phys_try_closure,summary_stat_try_closure,sim_param_closure)
  param_guess = make_vector_of_sim_param(sim_param_closure)
  evaluate_model_scalar_ret( param_guess)
end

end

