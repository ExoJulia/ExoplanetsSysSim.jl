# ExoplanetsSysSim/src/summary_statistics.jl
## (c) 2015 Eric B. Ford

using Statistics

mutable struct CatalogSummaryStatistics
  stat::Dict{String,Any}          # For storing summary statistics
  cache::Dict{String,Any}         # For caching data that's not a summary statistic
end

function CatalogSummaryStatistics()    # Constructor for empty CatalogSummaryStatistics objects for creating globals to be used by closure
  CatalogSummaryStatistics( Dict{String,Any}(), Dict{String,Any}() )
end

function calc_summary_stats_sim_pass_one_demo(cat_obs::KeplerObsCatalog, cat_phys::KeplerPhysicalCatalog, param::SimParam )      # Version for simulated data, since includes cat_phys
  ssd = Dict{String,Any}()
  cache = Dict{String,Any}()

  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1
  idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = findall(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target[idx_tranets] )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target[idx_tranets] )

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations

  cache["num_tranets"] = num_tranets
  cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we won't use these for computing a distance or probability
  #cache["idx_n_tranets"] = idx_n_tranets

  expected_num_detect = 0.0
  expected_num_sys_n_tranets = zeros(max_tranets_in_sys)
  for i in idx_tranets
    for j in 1:num_planets(cat_obs.target[i])
      p_tr_and_det = prob_detect(cat_obs.target[i].prob_detect,j)   # WARNING: Check why not using cat_phys here?
      expected_num_detect += p_tr_and_det
    end
    for k in 1:max_tranets_in_sys
      expected_num_sys_n_tranets[k] += prob_detect_n_planets(cat_obs.target[i].prob_detect,k)  # WARNING: Check why not use cat_phys here?
    end
  end
  ssd["expected planets detected"] = expected_num_detect
  ssd["num_sys_tranets"] = expected_num_sys_n_tranets
  ssd["num targets"] = get_int(param,"num_targets_sim_pass_one")
  #println("expected planets = ",expected_num_detect,", num_sys_tranets = ",expected_num_sys_n_tranets,", num targets = ",ssd["num targets"])

  # Arrays to store values for each tranet
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  weight_list = zeros(num_tranets)

  tr_id = 1   # tranet id
  for i in idx_tranets     # For each target with at least one tranet
     targ = cat_obs.target[i]
     for j in 1:min(length(targ.obs),max_tranets_in_sys)   # For each tranet around that target (but truncated if too many tranets in one system)
         #println("# i= ",i," j= ",j," tr_id= ",tr_id)
         period_list[tr_id] = targ.obs[j].period
         depth_list[tr_id] = targ.obs[j].depth
         # (s,p) = targ.phys_id[j]
         # ptr = calc_transit_prob_single(cat_phys.target[i],s,p)   # WARNING: Could access physical catalog, rather than observed catalog, but obviously that's dangerous for observations.
	       weight_list[tr_id] = prob_detect(cat_obs.target[i].prob_detect,j)
         tr_id += 1
      end
   end
  ssd["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  ssd["depth list"] = depth_list
  ssd["weight list"] = weight_list

  idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 && weight_list[i]>0.0 for i in 1:length(period_list) ]
  log_period_list = log10.(period_list[idx_good])
  log_depth_list = log10.(depth_list[idx_good])
  weight_list = weight_list[idx_good]
  weight_sum = sum(weight_list)
  ssd["mean log10 P"]  =  mean_log_P = sum( weight_list .* log_period_list) / weight_sum                           # TODO TEST: Check that these four weighted mean and stddevs are working properly
  ssd["std log10 P"]  =  sum( weight_list .* (log_period_list.-mean_log_P).^2 ) / weight_sum
  ssd["mean log10 depth"]  =  mean_log_depth = sum( weight_list .* log_depth_list) / weight_sum
  ssd["std log10 depth"]  =  sum( weight_list .* (log_depth_list.-mean_log_depth).^2 ) / weight_sum

  return CatalogSummaryStatistics(ssd, cache)
end


function calc_summary_stats_obs_demo(cat_obs::KeplerObsCatalog, param::SimParam )      # Version for observed data, thus no use of cat_phys
  ssd = Dict{String,Any}()
  cache = Dict{String,Any}()

  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")                  # Demo that simulation parameters can specify how to evalute models, too
  idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet
  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ [] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = findall(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target[idx_tranets] )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target[idx_tranets] )
  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  if ( length( findall(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target[idx_tranets] ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets = Int64(num_tranets)           # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations
  #println("# num_tranets= ",num_tranets)

  # QUERY:  Is there any reason to cache anything for the real observations?  We only need to do this once, so might as well use one pass to simplicity.
  #cache["num_tranets"] = num_tranets
  #cache["idx_tranets"] = idx_tranets                                   # We can save lists of indices to summary stats for pass 2, even though we won't use these for computing a distance or probability
  cache["idx_n_tranets"] = idx_n_tranets

  ssd["planets detected"] = num_tranets                                 # WARNING: Note that we'll comparing two different things for simulated and real data during pass 1 (expected planets detected)
  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  ssd["num_sys_tranets"] = num_sys_tranets
  ssd["num targets"] = get_int(param,"num_kepler_targets")

  # Arrays to store values for each tranet
  period_list = zeros(num_tranets)
  depth_list = zeros(num_tranets)
  weight_list = ones(num_tranets)

   i = 1   # tranet id
   for targ in cat_obs.target[idx_tranets]                        # For each target with at least one tranet
     for j in 1:min(length(targ.obs),max_tranets_in_sys)          # For each tranet around that target (but truncated if too many tranets in one system)
         #println("# i= ",i," j= ",j)
         period_list[i] = targ.obs[j].period
         depth_list[i] = targ.obs[j].depth
         #weight_list[i] = 1.0
         i = i+1
      end
   end

  ssd["P list"] = period_list                                     # We can store whole lists, e.g., if we want to compute K-S distances
  ssd["depth list"] = depth_list
  ssd["weight list"] = weight_list

  idx_good = Bool[ period_list[i]>0.0 && depth_list[i]>0.0 for i in 1:length(period_list) ]
  log_period_list = log10.(period_list[idx_good])
  log_depth_list = log10.(depth_list[idx_good])
  ssd["mean log10 P"]  =  mean_log_P = mean(log_period_list)
  ssd["mean log10 depth"]  =  mean_log_depth = mean(log_depth_list)
  ssd["std log10 P"]  =  stdm(log_period_list,mean_log_P)
  ssd["std log10 depth"]  =  stdm(log_depth_list,mean_log_depth)

  return CatalogSummaryStatistics(ssd, cache)
end


# Just returns summary statistics passed, but provides a demo/hook for computing more expensive summary statistics if a model is good enouguh to be worth the extra time.
function calc_summary_stats_sim_pass_two_demo(cat_obs::KeplerObsCatalog, cat_phys::KeplerPhysicalCatalog, ss::CatalogSummaryStatistics, param::SimParam )
  return ss
end

function test_summary_statistics(cat_obs::KeplerObsCatalog, cat_phys::KeplerPhysicalCatalog, sim_param::SimParam)
  ss = calc_summary_stats_sim_pass_one_demo(cat_obs,cat_phys,sim_param)
  #println("len (ss pass 1)= ",length(collect(keys(ss1.stat))))
  ss = calc_summary_stats_sim_pass_two_demo(cat_obs,cat_phys,ss,sim_param)
  #println("len (ss pass 1)= ",length(collect(keys(ss1.stat))), "...   len (ss pass 2)= ",length(collect(keys(ss2.stat))) )
  ss = calc_summary_stats_obs_demo(cat_obs,sim_param)   # So tests can compare to simulated observed catalog
  return ss
end
