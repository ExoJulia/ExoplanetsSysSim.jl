## ExoplanetsSysSim/examples/tess/abc_setup.jl
## (c) 2019 Eric B. Ford & Danley C. Hsu
# Collection of functions which specific ABC simulation parameters

module EvalSysSimModel
export setup, get_param_vector, get_ss_obs
export gen_data, calc_summary_stats, calc_distance, is_valid_uniform, is_valid_beta, is_valid_dirichlet, normalize_dirch
using ExoplanetsSysSim
using CSV
using DataFrames
using ApproximateBayesianComputing
const ABC = ApproximateBayesianComputing
include(joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")),"examples","tess", "tess_binrates_func.jl"))

sim_param_closure = SimParam()
summary_stat_ref_closure =  CatalogSummaryStatistics()

function is_valid_uniform(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    local rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    limitP::Array{Float64,1} = get_any(sim_param_closure, "p_lim_arr", Array{Float64,1})
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< sum(rate_tab, dims=1)')
        return false
    end
    return true
end

function is_valid_beta(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    local rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    limitP::Array{Float64,1} = get_any(sim_param_closure, "p_lim_arr", Array{Float64,1})
    local bin_size_factor::Float64 = get_real(sim_param_closure, "bin_size_factor")
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< (bin_size_factor*3*[log(limitP[i+1]/limitP[i])/log(2) for i in 1:length(limitP)-1].*sum(rate_tab, dims=1)'))
        return false
    end
    return true
end

function is_valid_dirichlet(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    local rate_tab::Array{Float64,2} = get_any(sim_param_closure, "obs_par", Array{Float64,2})
    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    #const lambda = sum_kbn(rate_tab)
    if any(x -> x < 0., rate_tab) || any([floor(3*log(limitP[i+1]/limitP[i])/log(2)) for i in 1:length(limitP)-1] .< rate_tab[1,:])
        return false
    end
    return true
end

function normalize_dirch(param_vector::Vector{Float64})
    global sim_param_closure
    local p_dim = length(get_any(sim_param_closure, "p_lim_arr", Array{Float64,1}))-1
    local r_dim = length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1

    for i in 1:p_dim
        param_vector[((i-1)*(r_dim+1)+2):((i-1)*(r_dim+1)+(r_dim+1))] ./= sum(param_vector[((i-1)*(r_dim+1)+2):((i-1)*(r_dim+1)+(r_dim+1))])
    end

    update_sim_param_from_vector!(param_vector,sim_param_closure)
    return param_vector
end

function gen_data(param_vector::Vector{Float64})
    global sim_param_closure
    update_sim_param_from_vector!(param_vector,sim_param_closure)
    cat_phys = generate_tess_physical_catalog(sim_param_closure)
    #cat_phys_cut = ExoplanetsSysSim.generate_obs_targets(cat_phys, sim_param_closure)
    #cat_obs = ExoplanetsSysSim.observe_tess_targets_single_obs(cat_phys_cut, sim_param_closure)
    cat_obs = ExoplanetsSysSim.observe_tess_targets_sky_avg(cat_phys, sim_param_closure)
    return cat_obs
end

# TODO OPT: Eventually, could adapt ABC.jl to use distance from first pass to decide if should compute additional summary statistics
function calc_summary_stats(cat::TESSObsCatalog)
    global sim_param_closure
    sum_stat = calc_summary_stats_obs_binned_rates(cat, sim_param_closure, obs_skyavg = true)
    return sum_stat
end

function calc_distance(sum_stat_obs::CatalogSummaryStatistics,sum_stat_sim::CatalogSummaryStatistics, n::Integer = 0)
    global sim_param_closure
    dist1 = calc_distance_vector_binned(sum_stat_obs,sum_stat_sim, 1, sim_param_closure)
    num_available = length(dist1)
    num_to_use = n>0 ? min(n,num_available) : num_available
    return calc_scalar_distance(dist1[1:num_to_use])
end

function read_star_table_tess()
    global sim_param_closure
    stlr = DataFrame()
    for i=1:27
        fname="tesstargets/TESS_targets_S" * lpad(i, 3, "0") * ".csv" # this is specific to my machine!
        stlr = vcat(stlr, CSV.read(fname))
    end
    add_param_active(sim_param_closure, "read_stellar_catalog", stlr)
    sim_param_closure
end

function setup(prior_choice::String, bin_size_factor::Float64)
    global sim_param_closure = setup_sim_param_tessbinrates()
    add_param_fixed(sim_param_closure,"bin_size_factor",bin_size_factor)
    if prior_choice == "dirichlet"
        sim_param_closure = set_test_param_total(sim_param_closure)
        add_param_fixed(sim_param_closure,"generate_num_planets",generate_num_planets_binrates_dirichlet)
        if (length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1) > 1
            add_param_fixed(sim_param_closure,"generate_period_and_sizes", generate_period_and_sizes_binrates_dirichlet)
        end
    elseif prior_choice == "beta"
        sim_param_closure = set_test_param(sim_param_closure)
        add_param_fixed(sim_param_closure,"generate_num_planets",generate_num_planets_binrates_beta)
        add_param_fixed(sim_param_closure,"generate_period_and_sizes", generate_period_and_sizes_binrates_beta)
    elseif prior_choice == "uniform"
        sim_param_closure = set_test_param(sim_param_closure)
    else
        println("# Invalid prior given!")
        quit()
    end
    
    ### Use simulated planet candidate catalog data
    # df_star = setup_star_table_christiansen(sim_param_closure)
    # println("# Finished reading in stellar data")
    # add_param_fixed(sim_param_closure,"num_tess_targets",1000000)  # For "observed" catalog
    # cat_obs = simulated_read_tess_observations(sim_param_closure)
    # println("# Finished setting up simulated true catalog")
    ###
    
    ### Use real planet candidate catalog data
    df_toi, usable_toi = read_toi_catalog(sim_param_closure)
    println("# Finished reading in TOI data")  
    df_star = setup_star_table_tic(sim_param_closure)
    println("# Finished reading in stellar data")
    cat_obs = setup_actual_pc_catalog_tess(df_star, df_toi, usable_toi, sim_param_closure)
    println("# Finished setting up true catalog")
    ###
    
    global summary_stat_ref_closure = calc_summary_stats_obs_binned_rates(cat_obs, sim_param_closure, trueobs_cat = true)
end

get_param_vector() = make_vector_of_sim_param(sim_param_closure)
get_ss_obs() = summary_stat_ref_closure

function set_simparam_ss(sim_param::ExoplanetsSysSim.SimParam, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics)
    global sim_param_closure = sim_param
    global summary_stat_ref_closure = ss_true
end

end  # module EvalSysSimModel

#include(joinpath(Pkg.dir("ABC"),"src/composite.jl"))

module SysSimABC
export setup_abc, run_abc, run_abc_largegen, setup_abc_p2
using Distributions, Random, Distributed
using ApproximateBayesianComputing
const ABC = ApproximateBayesianComputing
import ApproximateBayesianComputing.CompositeDistributions.CompositeDist
import ApproximateBayesianComputing.TransformedBetaDistributions.LinearTransformedBeta

#using Compat
import ExoplanetsSysSim
import ..EvalSysSimModel
include(joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")),"examples","tess", "tess_binrates_func.jl"))
include(joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")),"examples","tess", "beta_proposal.jl"))

function setup_abc(num_dist::Integer = 0; prior_choice::String = "uniform", bin_size_factor::Float64 = 1.5)
    EvalSysSimModel.setup(prior_choice, bin_size_factor)

    limitP::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "p_lim_arr", Array{Float64,1})
    limitR::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "r_lim_arr", Array{Float64,1})
    limitR_full::Array{Float64,1} = get_any(EvalSysSimModel.sim_param_closure, "r_lim_full", Array{Float64,1})
    local r_dim = length(limitR)-1
    
    prior_arr = ContinuousDistribution[]
    ss_obs_table = EvalSysSimModel.get_ss_obs().stat["planets table"]

    if prior_choice == "dirichlet"
        weights_arr = [log(limitR[j+1]/limitR[j]) for j in 1:r_dim]/minimum([log(limitR_full[k+1]/limitR_full[k]) for k in 1:(length(limitR_full)-1)])
        for i in 1:(length(limitP)-1)
            max_in_col = 3*log(limitP[i+1]/limitP[i])/log(2)
            lambda_col = Uniform(0.0, max_in_col)
            prior_arr = vcat(prior_arr, lambda_col)
            if r_dim > 1
                dirch_dist = Dirichlet(weights_arr)
                prior_arr = vcat(prior_arr, dirch_dist)
            end
        end
    elseif prior_choice == "beta"
        for i in 1:(length(limitP)-1)
            for j in 1:r_dim
                r_ind = findfirst(x -> x == limitR[j], limitR_full)
                beta_dist = Beta(log(limitR[j+1]/limitR[j]), sum(log.([getindex(limitR_full, x) for x = 2:length(limitR_full) if x != r_ind+1] ./ [getindex(limitR_full, x) for x = 1:length(limitR_full)-1 if x != r_ind])))
                prior_arr = vcat(prior_arr, beta_dist)
            end
        end
    else
        for i in 1:(length(limitP)-1)
            max_in_col = bin_size_factor*log(limitP[i+1]/limitP[i])/log(2)
            for j in 1:r_dim
                uniform_dist = Uniform(0.0, max_in_col*log(limitR[j+1]/limitR[j])/log(2))
                prior_arr = vcat(prior_arr, uniform_dist)
            end
        end
    end

    param_prior = CompositeDist(prior_arr)
    in_parallel = nworkers() > 1 ? true : false

    calc_distance_ltd(sum_stat_obs::ExoplanetsSysSim.CatalogSummaryStatistics,sum_stat_sim::ExoplanetsSysSim.CatalogSummaryStatistics) = EvalSysSimModel.calc_distance(sum_stat_obs,sum_stat_sim,num_dist)

    global abc_plan = ABC.abc_pmc_plan_type(EvalSysSimModel.gen_data,EvalSysSimModel.calc_summary_stats, calc_distance_ltd, param_prior, make_proposal_dist=make_proposal_dist_multidim_beta, is_valid=EvalSysSimModel.is_valid_uniform, num_part=200, num_max_attempt=400, num_max_times=30, epsilon_init=9.9e99, target_epsilon=1.0e-100, in_parallel=in_parallel, adaptive_quantiles = false, epsilon_reduction_factor=0.9, tau_factor=1.5);

    if prior_choice == "dirichlet" && r_dim > 1
        abc_plan.make_proposal_dist = make_proposal_dist_multidim_beta_dirichlet
        abc_plan.is_valid = EvalSysSimModel.is_valid_dirichlet
        abc_plan.normalize = EvalSysSimModel.normalize_dirch
    elseif prior_choice == "beta"
        abc_plan.is_valid = EvalSysSimModel.is_valid_beta
    end
    
    return abc_plan
end

function setup_abc_p2(abc_plan::ABC.abc_pmc_plan_type)
    abc_plan.tau_factor = 2.0
    abc_plan.num_max_times = 200
    ExoplanetsSysSim.add_param_fixed(EvalSysSimModel.sim_param_closure,"num_targets_sim_pass_one",10000)
    return abc_plan
end
    

function run_abc_largegen(abc_plan::ABC.abc_pmc_plan_type, pop::ABC.abc_population_type, ss_true::ExoplanetsSysSim.CatalogSummaryStatistics, epshist_targ::Float64; npart::Integer = 1000, num_dist::Integer = 0)

    abc_plan.num_max_times = 1

    println("# run_abc_largegen: ",EvalSysSimModel.sim_param_closure)
    sampler_largegen = abc_plan.make_proposal_dist(pop, abc_plan.tau_factor)
    theta_largegen = Array{Float64}(undef, size(pop.theta, 1), npart)
    weight_largegen = Array{Float64}(undef, npart)
    for i in 1:npart
        theta_val, dist_largegen, attempts_largegen = ABC.generate_theta(abc_plan, sampler_largegen, ss_true, epshist_targ)
        theta_largegen[:,i] = theta_val  
        prior_logpdf = Distributions.logpdf(abc_plan.prior,theta_val)
        sampler_logpdf = Distributions.logpdf(sampler_largegen, theta_val)
        weight_largegen[i] = exp(prior_logpdf-sampler_logpdf)
    end
    return theta_largegen, weight_largegen
end

function run_abc(abc_plan::ABC.abc_pmc_plan_type)
    #global sim_param_closure
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    #println("True catalog SS: ", ss_true)
    pop_out = ABC.run_abc(abc_plan,ss_true;verbose=true)
end

function run_abc(abc_plan::ABC.abc_pmc_plan_type, pop::ABC.abc_population_type)
    #global sim_param_closure
    dist_threshold = maximum(pop.dist)
    EvalSysSimModel.add_param_fixed(EvalSysSimModel.sim_param_closure,"minimum ABC dist skip pass 2",dist_threshold)
    println("# run_abc: ",EvalSysSimModel.sim_param_closure)
    ss_true = EvalSysSimModel.get_ss_obs()
    pop_out = ABC.run_abc(abc_plan,ss_true,pop;verbose=true)
end

end # module SysSimABC
