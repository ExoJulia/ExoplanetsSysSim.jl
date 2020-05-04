## ExoplanetsSysSim/examples/hsu_etal_2018/invdet_calc.jl
## (c) 2018 Danley C. Hsu
# Script for producing DR25 FGK planet candidate occurrence rate estimates
#    using both the inverse detection efficiency and the simple
#    Bayesian methods

using ExoplanetsSysSim
include(joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")),"examples","dr25_gaia_fgk", "dr25_binrates_func.jl"))

global sim_param_closure = setup_sim_param_dr25binrates()
sim_param_closure = set_test_param(sim_param_closure)

df_koi,usable_koi = read_koi_catalog(sim_param_closure)
println("# Finished reading in KOI data")
df_star = setup_star_table_dr25(sim_param_closure)
println("# Finished reading in stellar data") 
cat_obs = setup_actual_planet_candidate_catalog(df_star, df_koi, usable_koi, sim_param_closure)

#@time inv_det_simp_bayes(cat_obs, sim_param_closure)
@time simp_bayes(cat_obs, sim_param_closure)
