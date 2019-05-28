## ExoplanetsSysSims/src/abc_distance.jl
## (c) 2015 Eric B. Ford

dist_L1_fractional(x::Real, y::Real) =  (x==y==0.0) ? 0.0 : abs( (x-y)*2/(x+y) )
dist_L1_abs(x::Real, y::Real) =  abs(x-y)
dist_L2_fractional(x::Real, y::Real) =  (x==y==0.0) ? 0.0 : abs( (x-y)^2*2/(x^2+y^2) )
dist_L2_abs(x::Real, y::Real) =  (x-y)^2

# Library of functions to convert distance vector into a scalar distance
calc_scalar_distance_sum(d::Vector{Float64}) = sum(d)
calc_scalar_distance_rms(d::Vector{Float64}) = sqrt(sumsq(d))
calc_scalar_distance_max(d::Vector{Float64}) = maximum(d)

# How to combine two distances based on subsets of the full distance vector (should always be greater than either for two pass algorithm to work)
combine_scalar_distances_sum(d1::Float64, d2::Float64) = d1+d2
combine_scalar_distances_rms(d1::Float64, d2::Float64) = sqrt(d1*d1+d2*d2)
combine_scalar_distances_max(d1::Float64, d2::Float64) = max(d1,d2)

# Pick which one to be used
calc_scalar_distance(d::Vector{Float64}) = calc_scalar_distance_sum(d)
combine_scalar_distances(d1::Float64, d2::Float64) = combine_scalar_distances_sum(d1,d2)

# compute supremum of differences between empirical cdfs.
# Borrowed from JuliaStats/HypothesisTests.jl
function ksstats(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    n_x, n_y = length(x), length(y)
    sort_idx = sortperm([x; y])
    pdf_diffs = [ones(n_x)/n_x; -ones(n_y)/n_y][sort_idx]
    cdf_diffs = cumsum(pdf_diffs)
    deltap = maximum(cdf_diffs)
    deltan = -minimum(cdf_diffs)
    delta = max(deltap,deltan)
    (n_x, n_y, deltap, deltan, delta)
end
# weighted version   # WARNING:  Function is untested
function ksstats(x::AbstractVector{T}, y::AbstractVector{S}, wx::AbstractVector{T}, wy::AbstractVector{T}) where {T<:Real, S<:Real}
    n_x, n_y = length(x), length(y)
    wx .*= 1.0/sum(wx)
    wy .*= 1.0/sum(wy)
    sort_idx = sortperm([x; y])
    pdf_diffs = [wx; -wy][sort_idx]
    cdf_diffs = cumsum(pdf_diffs)
    deltap = maximum(cdf_diffs)
    deltan = -minimum(cdf_diffs)
    delta = max(deltap,deltan)
    (n_x, n_y, deltap, deltan, delta)  # should the first two values returned here be sum(wx) and sum(wy) before normalizing?  For now, we ignore all but the final delta anyway.
end

dist_KS(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}  = ksstats(x,y)[5]
dist_KS(x::AbstractVector{T}, y::AbstractVector{S}, wx::AbstractVector{T}, wy::AbstractVector{T}) where {T<:Real, S<:Real} = ksstats(x,y,wx,wy)[5]

# lambda:  rate for Poisson process, i.e., expected value for number of events
# k:  number of events observed

function calc_num_events_to_maximize_poisson_pdf(lambda::Real)
  kstar = floor(Int64,lambda)
  if lambda > kstar+1
     kstar +=1
  end
  return kstar
end

function log_prob_poisson_num_events_given_rate(lambda::Real,k::Real)
   k*log(lambda)-lambda-lgamma(k+1)
end

function delta_log_prob_poisson_num_events_given_rate(lambda::Real, k::Integer, kstar::Integer = calc_num_events_to_maximize_poisson_pdf(lambda) )
  delta = (k-kstar)*log(lambda)
  if kstar>=k+1
	delta += sum(log.((k+1):kstar))
  elseif k>=kstar+1
     delta -= sum(log.((kstar+1):k))
  end
  return delta
end

function distance_poisson_likelihood(lambda::Real, k::Integer)
  kstar = calc_num_events_to_maximize_poisson_pdf(lambda)
  khi = round(Int64,lambda+sqrt(lambda))
  delta_logprob_one_sigma_hi = delta_log_prob_poisson_num_events_given_rate(lambda,khi,kstar)
  delta_logp_k = delta_log_prob_poisson_num_events_given_rate(lambda,k,kstar)
  delta_logp_k/delta_logprob_one_sigma_hi
end

function distance_poisson_draw(lambda::Real, k::Integer)
  d =  Distributions.Poisson(lambda)
  simulated_number_of_detections = rand(d)
  abs( simulated_number_of_detections -k)
end


function distance_sum_of_bernoulli_draws(num_pl_obs::Integer, num_targets_obs::Integer, prob_detect_list::Vector{TReal}, num_targets_sim::Integer, bin_match_list::Vector{TInt}) where {TReal<:Real, TInt<:Integer}
   @assert(0<=num_pl_obs<=num_targets_obs)
   num_pl_match = length(bin_match_list)
   @assert(0<=length(bin_match_list))

   num_detect_sim = 0
   if num_pl_match >= 1
      num_draws_all = min(max(1,floor(Int64, num_targets_obs/num_targets_sim)),1000)
      @assert(1<=num_draws_all<=1000)
      for i in 1:num_pl_match
         pl_id = bin_match_list[i]
         @assert(1<=pl_id<=length(prob_detect_list))
         prob_detect = min(prob_detect_list[pl_id],1.0)
         num_detect_sim += sum(rand(Bernoulli(prob_detect),num_draws_all))
      end
      # If number of targets observed is not a multiple of number of targets simulated, then pick a random set to make total number of draws equal (as long as there are some planets to choose from)
      for i in (num_pl_match*num_draws_all+1):(floor(Int64,num_pl_match*num_targets_obs/num_targets_sim))
          pl_id = bin_match_list[rand(1:num_pl_match)]
          prob_detect = min(prob_detect_list[pl_id],1.0)
          num_detect_sim += rand(Bernoulli(prob_detect))
      end
   end
   distance = dist_L2_abs(num_pl_obs/num_targets_obs, num_detect_sim/num_targets_obs)
   return distance, num_detect_sim
end

# compute Canberra distance.
function distance_canberra_modified(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    @assert length(x) == length(y)
    dist_sum = 0.0
    for i in 1:length(x)
        numer = abs(x[i]-y[i])
        denom = sqrt(x[i] + y[i])
        if denom == 0.0
            continue
        else
            dist_sum += numer/denom
        end
    end
    return dist_sum
end

# compute originalCanberra distance.
function distance_canberra_orig(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    @assert length(x) == length(y)
    dist_sum = 0.0
    for i in 1:length(x)
        numer = abs(x[i]-y[i])
        denom = abs(x[i]) + abs(y[i])
        if denom == 0.0
            continue
        else
            dist_sum += numer/denom
        end
    end
    return dist_sum
end

distance_canberra(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real} = distance_canberra_modified(x,y)  # WARNING: Defaults to modified Canberra distance.  At some point, update Danley's code to call distance_canberra_modified, so can remove this

# compute Cosine distance.
function distance_cosine(x::AbstractVector{T}, y::AbstractVector{S}) where {T<:Real, S<:Real}
    @assert length(x) == length(y)
    numer = 0.0
    denom_1 = 0.0
    denom_2 = 0.0
    for i in 1:length(x)
        numer += x[i]*y[i]
        denom_1 += x[i]^2
        denom_2 += y[i]^2
    end
    return numer/(sqrt(denom_1)*sqrt(denom_2))
end

function wasserstein_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2}, p::Integer) where { T1<:Real, T2<:Real}
  # code adapted from SciPy.stats._cdf_distance
  u_sorter = issorted(u_values) ? (1:length(u_values)) : sortperm(u_values)
  v_sorter = issorted(v_values) ? (1:length(v_values)) : sortperm(v_values)
  all_values = sort(vcat(u_values,v_values),alg=MergeSort)
  deltas = all_values[2:end] .- all_values[1:end-1]
  u_cdf = map(x->searchsortedlast(u_values[u_sorter],x)/ length(u_values),all_values[1:end-1])  
  v_cdf = map(x->searchsortedlast(v_values[v_sorter],x)/ length(v_values),all_values[1:end-1])
  if p==1
     return sum(abs.(u_cdf.-v_cdf).*deltas)
    elseif p==2
      return sqrt(sum((u_cdf.-v_cdf).^2 .* deltas))
    else
      return pow(sum(abs.(u_cdf.-v_cdf).^p .* deltas),1/p)
    end
end

earth_mover_distance(u_values::AbstractVector{T1}, v_values::AbstractVector{T2}) where { T1<:Real, T2<:Real} = wasserstein_distance(u_values,v_values,1)

# TODO USER SCI: IMPORTANT:  Replace the distance function with something well thought out for your particular scientific application.  See examples
function calc_distance_vector_demo(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, pass::Int64, sim_param::SimParam ; verbose::Bool = false)
  d = Array{Float64}(undef,0)
  if pass == 1
    if verbose
      println("# Summary 1, pass 1: ",summary1)
      println("# Summary 2, pass 1: ",summary2)
    end
    d = zeros(3)
    # Since observed and simulated catalogs can have different summary statistics for the number of planets, prefer detections if avaliable (e.g., after pass2), otherwise use expected (e.g., from pass 1)
    np1 = haskey(summary1.stat,"planets detected") ? summary1.stat["planets detected"] : summary1.stat["expected planets detected"]
    np2 = haskey(summary2.stat,"planets detected") ? summary2.stat["planets detected"] : summary2.stat["expected planets detected"]
    d[1] = dist_L1_abs(np1/summary1.stat["num targets"],np2/summary2.stat["num targets"])    #  Normalize so different statistics weighted appropriately and not dominated by this one
    #println("np1 = ",np1,", np2 = ",np2)
    #println("np1 (normalized) = ",np1/summary1.stat["num targets"],", np2 (normalized) = ",np2/summary2.stat["num targets"],", d[1] = ",d[1])
    #d[2] = dist_L1_abs(summary1.stat["mean log10 P"],summary2.stat["mean log10 P"])
    #d[3] = dist_L1_abs(summary1.stat["mean log10 depth"],summary2.stat["mean log10 depth"])
    #d[4] = dist_L1_abs(summary1.stat["std log10 P"],summary2.stat["std log10 P"])
    #d[5] = dist_L1_abs(summary1.stat["std log10 depth"],summary2.stat["std log10 depth"])
    #d[2] = dist_KS(summary1.stat["P list"], summary2.stat["P list"])
    #d[3] = dist_KS(summary1.stat["depth list"], summary2.stat["depth list"])
    # EDITS for christiansen-single-bin
    #d[2] = dist_KS(summary1.stat["P list"], summary2.stat["P list"],summary1.stat["weight list"],summary2.stat["weight list"])
    #d[3] = dist_KS(summary1.stat["depth list"], summary2.stat["depth list"],summary1.stat["weight list"],summary2.stat["weight list"])
    # END EDITS
  elseif pass == 2
    max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
    d = zeros(max_tranets_in_sys)
    for n in 1:max_tranets_in_sys
      d[n] = n*dist_L1_abs(summary1.stat["num_sys_tranets"][n]/summary1.stat["num targets"],summary2.stat["num_sys_tranets"][n]/summary2.stat["num targets"])
    end
  else
    println("# calc_distance_vector_demo doesn't know what to do for pass= ", pass)
  end
  return d
end

function test_abc_distance(cat_obs::KeplerObsCatalog, cat_phys::KeplerPhysicalCatalog, sim_param::SimParam)
  ss_pass1 = calc_summary_stats_sim_pass_one_demo(cat_obs,cat_phys,sim_param)
  ss_pass2 = calc_summary_stats_sim_pass_two_demo(cat_obs,cat_phys,ss_pass1,sim_param)
  d1 = calc_distance_vector_demo(ss_pass1,ss_pass1, 1, sim_param)
  d2 = calc_distance_vector_demo(ss_pass2,ss_pass2, 2, sim_param)
end
