## ExoplanetsSysSim/src/ExoplanetsSysSim.jl
## (c) 2015 Eric B. Ford

using DataFrames # needed outside module, so JLD/JLD2 load has the type

module ExoplanetsSysSim

# Packages to be used                 # Note: Tried to keep to a minimum in core package to help with maintainability
import Compat: @compat #, readstring, is_windows
if VERSION < v"0.7"
  import Compat: UTF8String, ASCIIString
end
if VERSION >= v"0.5-"
  using Combinatorics
end
using Distributions
using PDMats
# using StatsFuns
using DataFrames
if VERSION < VersionNumber(0,3,20)
  using Docile
end

using CORBITS

# Includes & associated  exports for types, then generic functions, then demo functions
include("constants.jl")
#require(joinpath(Pkg.dir("ExoplanetsSysSim"), "src", "constants.jl"))
export SimParam
export SimulationParameters
export add_param_fixed, add_param_active, set_active, set_inactive, is_active, update_param, get, get_real, get_int, get_function, get_any, haskey, make_vector_of_sim_param, update_sim_param_from_vector!, make_vector_of_active_param_keys, get_range_for_sim_param
export setup_sim_param_demo
include("simulation_parameters.jl")
using ExoplanetsSysSim.SimulationParameters
export Orbit
include("orbit.jl")
export Planet
include("planet.jl")
export LimbDarkeningParamAbstract, LimbDarkeningParamLinear, LimbDarkeningParamQuadratic, LimbDarkeningParam4thOrder
include("limb_darkening.jl")
export depth_at_midpoint, ratio_from_depth
export StarAbstract, Star, SingleStar, BinaryStar
include("star.jl")
export flux, mass
export generate_stars
export PlanetarySystemAbstract, PlanetarySystemSingleStar, PlanetarySystem
include("planetary_system.jl")
#include("corbits.jl")
#export prob_of_transits_approx
include("window_function.jl")
export WindowFunction
export setup_window_function, get_window_function_data, get_window_function_id, eval_window_function, setup_OSD_interp#, cdpp_vs_osd
include("stellar_table.jl")
export  StellarTable
export setup_star_table, star_table, num_usable_in_star_table, set_star_table, star_table_has_key
export KeplerTarget
export num_planets, generate_kepler_target_from_table, generate_kepler_target_simple
include("target.jl")
export KeplerTargetObs
include("transit_observations.jl")
export semimajor_axis
#export setup_koi_table, koi_table, num_koi
export KeplerPhysicalCatalog, KeplerObsCatalog
include("kepler_catalog.jl")
export generate_kepler_physical_catalog, observe_kepler_targets_single_obs, observe_kepler_targets_sky_avg, simulated_read_kepler_observations, setup_actual_planet_candidate_catalog, read_koi_catalog
export CatalogSummaryStatistics
include("summary_statistics.jl")
export calc_summary_stats_obs_demo, calc_summary_stats_sim_pass_one_demo, calc_summary_stats_sim_pass_two_demo
include("abc_distance.jl")
export dist_L1_fractional, dist_L1_abs, dist_L2_fractional, dist_L2_abs, calc_scalar_distance, combine_scalar_distances, distance_poisson_likelihood, distance_poisson_draw, distance_sum_of_bernoulli_draws, distance_canberra, distance_cosine
export calc_distance_vector_demo
export TestEvalModel
include("eval_model.jl")   # Also includes macros to help write eval model using different variables for closures
export test_eval_model
using ExoplanetsSysSim.TestEvalModel
export SysSimIO
include("io.jl")
using ExoplanetsSysSim.SysSimIO
export save_sim_param, save_sim_results, load_sim_param, load_distances, load_summary_stats

end # module
