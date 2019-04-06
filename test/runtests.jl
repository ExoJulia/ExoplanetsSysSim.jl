using Test

# write your own tests here
@test 1 == 1

using DataFrames
using ExoplanetsSysSim

function run_constructor_tests()
  ExoplanetsSysSim.SimulationParameters.test_sim_param_constructors()
  sim_param = ExoplanetsSysSim.setup_sim_param_demo()
  ExoplanetsSysSim.test_orbit_constructors()
  ExoplanetsSysSim.test_planet_constructors(sim_param)
  ExoplanetsSysSim.test_star_constructors(sim_param)
  ExoplanetsSysSim.test_planetary_system_constructors(sim_param)
  ExoplanetsSysSim.test_corbits()
  ExoplanetsSysSim.test_target(sim_param)
  ExoplanetsSysSim.test_transit_observations(sim_param)
  (cat_phys, cat_obs) = ExoplanetsSysSim.test_catalog_constructors(sim_param)
  ExoplanetsSysSim.test_summary_statistics(cat_obs, cat_phys, sim_param)
  ExoplanetsSysSim.test_abc_distance(cat_obs, cat_phys, sim_param)
  return 0
end

@test run_constructor_tests() == 0  # Just tests that the basic elements compile and run  # TODO: Write tests that will be useful in diagnosing any bugs
