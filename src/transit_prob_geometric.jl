## ExoplanetsSysSim/src/transit_prob_geometric.jl
## (c) 2015 Eric B. Ford

# Transit probability expressions if treat planets separately
function calc_transit_prob_single_planet_approx(P::Float64, Rstar::Float64, Mstar::Float64)
  return min(rsol_in_au*Rstar/semimajor_axis(P,Mstar), 1.0)
end

#function calc_transit_prob_single_planet_obs_ave(ps::PlanetarySystemAbstract, pl::Integer)
function calc_transit_prob_single_planet_obs_ave(ps::PlanetarySystem{StarT}, pl::Integer) where {StarT<:StarAbstract}

  ecc::Float64 = ps.orbit[pl].ecc
 a::Float64 = semimajor_axis(ps,pl)
 Rstar::Float64 = rsol_in_au*ps.star.radius
 return min(Rstar/(a*(1-ecc)*(1+ecc)), 1.0)
end
calc_transit_prob_single_planet_obs_ave(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_prob_single_planet_obs_ave(t.sys[s], p)

#=
# WARNING: This knows about e and w, but still returns a fraction rather than a 0 or 1.  Commented out for now, so no one uses it accidentally until we figure out why it was this way
function calc_transit_prob_single_planet_one_obs(ps::PlanetarySystemAbstract, pl::Integer)
 where {StarT<:StarAbstract}
  ecc::Float64 = ps.orbit[pl].ecc
 a::Float64 = semimajor_axis(ps,pl)
 Rstar::Float64 = rsol_in_au*ps.star.radius
 return min(Rstar*(1+ecc*sin(ps.orbit[pl].omega))/(a*(1-ecc)*(1+ecc)), 1.0)
end
calc_transit_prob_single_planet_one_obs(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_prob_single_planet_one_obs(t.sys[s], p)
=#

# WARNING: Assumes that planets with b>1 won't be detected/pass vetting
#function does_planet_transit(ps::PlanetarySystemAbstract, pl::Integer)
function does_planet_transit(ps::PlanetarySystem{StarT}, pl::Integer) where {StarT<:StarAbstract}
    ecc::Float64 = ps.orbit[pl].ecc
   incl::Float64 = ps.orbit[pl].incl
   a::Float64 = semimajor_axis(ps,pl)
   Rstar::Float64 = rsol_in_au*ps.star.radius
   if (Rstar >= (a*(1-ecc)*(1+ecc))/(1+ecc*sin(ps.orbit[pl].omega))*cos(incl))
     return true
   else
     return false
   end
end

#function corbits_placeholder_obs_ave( ps::PlanetarySystemSingleStar, use_pl::Vector{Cint} )    # Might be useful for someone to test w/o CORBITS
function corbits_placeholder_obs_ave( ps::PlanetarySystem{StarT}, use_pl::Vector{Cint} ) where {StarT<:StarAbstract}    # Might be useful for someone to test w/o CORBITS
  n = num_planets(ps)
  prob = 1.0
  for p in 1:n
     ptr = calc_transit_prob_single_planet_obs_ave(ps,p)
     prob *= (use_pl[p]==1) ? ptr : 1.0-ptr
     #=
     if(use_pl[p]==1)
       prob *= calc_transit_prob_single_planet_obs_ave(ps,p)
     else
       prob *= 1.0-calc_transit_prob_single_planet_obs_ave(ps,p)
     end
     =#
  end
  return prob
end

#function calc_impact_parameter(ps::PlanetarySystemSingleStar, pl::Integer)
function calc_impact_parameter(ps::PlanetarySystem{StarT}, pl::Integer) where {StarT<:StarAbstract}
      one_minus_e2 = (1-ps.orbit[pl].ecc)*(1+ps.orbit[pl].ecc)
      a_semimajor_axis = semimajor_axis(ps,pl)
      b = a_semimajor_axis *cos(ps.orbit[pl].incl)/(ps.star.radius*rsol_in_au)
      b *= one_minus_e2/(1+ps.orbit[pl].ecc*sin(ps.orbit[pl].omega))
      b = abs(b)
end

#function prob_combo_transits_one_obs( ps::PlanetarySystemSingleStar, use_pl::Vector{Cint} )
function prob_combo_transits_one_obs( ps::PlanetarySystem{StarT}, use_pl::Vector{Cint} )    where {StarT<:StarAbstract}
  n = num_planets(ps)
  for p in 1:n
      #one_minus_e2 = (1-ps.orbit[p].ecc)*(1+ps.orbit[p].ecc)
      #a_semimajor_axis = semimajor_axis(ps,p)
      #b = a_semimajor_axis *cos(ps.orbit[p].incl)/ps.star.radius
      #b *= one_minus_e2/(1+ps.orbit[p].ecc*sin(ps.orbit[p].omega))
      b = calc_impact_parameter(ps, p)
      if ! ( (b<=1.0 && use_pl[p]==1) || (b> 1.0 && use_pl[p]!=1) )
        return 0.0
      end
  end
  return 1.0
end

struct prob_combo_transits_obs_ave_workspace_type
  a::Vector{Cdouble}
  r::Vector{Cdouble}
  ecc::Vector{Cdouble}
  Omega::Vector{Cdouble}
  omega::Vector{Cdouble}
  inc::Vector{Cdouble}
  function prob_combo_transits_obs_ave_workspace_type(n::Integer)
     @assert(1<=n<=100)
     new( Array{Cdouble}(n), Array{Cdouble}(n), Array{Cdouble}(n), Array{Cdouble}(n), Array{Cdouble}(n), Array{Cdouble}(n) )
  end
end

#=
# Attempt to reduce memory allocations.  It does that, but no noticable speed improvement, so I'm leaving it commented out for now.
const global corbits_max_num_planets_per_system = 20
prob_combo_transits_obs_ave_workspace =  prob_combo_transits_obs_ave_workspace_type(corbits_max_num_planets_per_system)

#function prob_combo_transits_obs_ave( ps::PlanetarySystemSingleStar, use_pl::Vector{Cint}; print_orbit::Bool = false)
function prob_combo_transits_obs_ave( ps::PlanetarySystem{StarT}, use_pl::Vector{Cint}; print_orbit::Bool = false) where {StarT<:StarAbstract}
  n = num_planets(ps)
  @assert(n<=corbits_max_num_planets_per_system)
  for i in 1:n
       prob_combo_transits_obs_ave_workspace.a[i] = semimajor_axis(ps,i)
       prob_combo_transits_obs_ave_workspace.r[i] = ps.planet[i].radius * rsol_in_au
       prob_combo_transits_obs_ave_workspace.ecc[i] = ps.orbit[i].ecc
       prob_combo_transits_obs_ave_workspace.Omega[i] = ps.orbit[i].asc_node
       prob_combo_transits_obs_ave_workspace.omega[i] =ps.orbit[i].omega
       prob_combo_transits_obs_ave_workspace.inc[i] = ps.orbit[i].incl
  end

  r_star = convert(Cdouble,ps.star.radius *  rsol_in_au )
  prob = prob_of_transits_approx(prob_combo_transits_obs_ave_workspace.a, r_star, prob_combo_transits_obs_ave_workspace.r, prob_combo_transits_obs_ave_workspace.ecc, prob_combo_transits_obs_ave_workspace.Omega, prob_combo_transits_obs_ave_workspace.omega, prob_combo_transits_obs_ave_workspace.inc, use_pl)
  #prob = prob_of_transits_approx(a, r_star, r, ecc, Omega, omega, inc, use_pl)

  if print_orbit
  println("# a = ", prob_combo_transits_obs_ave_workspace.a)
  println("# r_star = ", r_star)
  println("# r = ", prob_combo_transits_obs_ave_workspace.r)
  println("# ecc = ", prob_combo_transits_obs_ave_workspace.ecc)
  println("# Omega = ", prob_combo_transits_obs_ave_workspace.Omega)
  println("# omega = ", prob_combo_transits_obs_ave_workspace.omega)
  println("# inc = ", prob_combo_transits_obs_ave_workspace.inc)
  println("# use_pl = ", use_pl)
  println("")
  end
  return prob
end
=#

#function prob_combo_transits_obs_ave( ps::PlanetarySystemSingleStar, use_pl::Vector{Cint}; print_orbit::Bool = false)
function prob_combo_transits_obs_ave( ps::PlanetarySystem{StarT}, use_pl::Vector{Cint}; print_orbit::Bool = false) where {StarT<:StarAbstract}
  n = num_planets(ps)
  a =  Cdouble[ semimajor_axis(ps,i) for i in 1:n ]
  r_star = convert(Cdouble,ps.star.radius *  rsol_in_au )
  r = Cdouble[ ps.planet[i].radius * rsol_in_au for i in 1:n ]
  ecc = Cdouble[ ps.orbit[i].ecc for i in 1:n ]
  Omega = Cdouble[ ps.orbit[i].asc_node for i in 1:n ]
  omega = Cdouble[ ps.orbit[i].omega for i in 1:n ]
  inc = Cdouble[ ps.orbit[i].incl for i in 1:n ]
  #use_pl = Cint[0 for i in 1:n]
  #for i in 1:length(combo)
  #  use_pl[i] = 1
  #end
  prob = prob_of_transits_approx(a, r_star, r, ecc, Omega, omega, inc, use_pl)

  if print_orbit
  println("# a = ", a)
  println("# r_star = ", r_star)
  println("# r = ", r)
  println("# ecc = ", ecc)
  println("# Omega = ", Omega)
  println("# omega = ", omega)
  println("# inc = ", inc)
  println("# use_pl = ", use_pl)
  println("")
  end
  return prob
end



@compat abstract type SystemDetectionProbsAbstract  end
@compat abstract type SystemDetectionProbsTrait end
@compat abstract type SkyAveraged <: SystemDetectionProbsTrait end
@compat abstract type OneObserver <: SystemDetectionProbsTrait end
# Derived types will allow us to specialize depending on whether using sky-averaged values, values for actual geometry (both of which require the physical catalog), or estimates based on observed data

mutable struct SimulatedSystemDetectionProbs{T<:SystemDetectionProbsTrait} <: SystemDetectionProbsAbstract         # To be used for simulated systems where we can calculat everything
  # Inputs to CORBITS
  detect_planet_if_transits::Vector{Float64}       # Probability of detecting each planet, averaged over all observers for each planet individually, assumes b~U[0,1);  To be used in pass 1

  # Outputs from CORBITS
  pairwise::Matrix{Float64}                 # detection probability (incl geometry & detection probability) for each planet (diagonal) and each planet pair (off diagonal)
                                            # TODO OPT: Make matrix symmetric to save memory?
  n_planets::Vector{Float64}                # fraction of time would detect n planets (incl. geometry & detection probability)

  combo_detected::Vector{Vector{Int64}}     # List of combinations of planets drawn from full joint multi-transit probability
end
# Removed since not really need and previous had used alised typename and function
# Specialize so know whether the values are sky averaged or not
#typealias SkyAveragedSystemDetectionProbs SimulatedSystemDetectionProbs{SkyAveraged}
#typealias OneObserverSystemDetectionProbs SimulatedSystemDetectionProbs{OneObserver}
#SkyAveragedSystemDetectionProbs = SimulatedSystemDetectionProbs{SkyAveraged}
#OneObserverSystemDetectionProbs = SimulatedSystemDetectionProbs{OneObserver}

function SimulatedSystemDetectionProbs(traits::Type, p::Vector{Float64}; num_samples::Integer = 1)
  SimulatedSystemDetectionProbs{traits}( p, zeros(length(p),length(p)), zeros(length(p)),  fill(Array{Int64}(undef,0), num_samples) )
end

function SimulatedSystemDetectionProbs(traits::Type, n::Integer; num_samples::Integer = 1)
  SimulatedSystemDetectionProbs{traits}( ones(n), zeros(n,n), zeros(n), fill(Array{Int64}(undef,0), num_samples) )
end

SkyAveragedSystemDetectionProbs(p::Vector{Float64}; num_samples::Integer = 1) = SimulatedSystemDetectionProbs( SkyAveraged, p, num_samples=num_samples)
SkyAveragedSystemDetectionProbs(n::Integer; num_samples::Integer = 1) = SimulatedSystemDetectionProbs(SkyAveraged, n, num_samples=num_samples)
SkyAveragedSystemDetectionProbsEmpty() = SimulatedSystemDetectionProbs(SkyAveraged,0)

OneObserverSystemDetectionProbs(p::Vector{Float64}; num_samples::Integer = 1) = SimulatedSystemDetectionProbs( OneObserver, p, num_samples=num_samples)

OneObserverSystemDetectionProbs(n::Integer; num_samples::Integer = 1) = SimulatedSystemDetectionProbs(OneObserver, n, num_samples=num_samples)
OneObserverSystemDetectionProbsEmpty() = SimulatedSystemDetectionProbs(OneObserver,0)

# Functions common to various types of SystemDetectionProbs
num_planets(prob::SimulatedSystemDetectionProbs{T}) where T<:SystemDetectionProbsTrait = length(prob.detect_planet_if_transits)
function prob_detect_if_transits(prob::SimulatedSystemDetectionProbs{T}, pl_id::Integer) where T<:SystemDetectionProbsTrait
  @assert 1<=pl_id<=length(prob.detect_planet_if_transits)
  prob.detect_planet_if_transits[pl_id]
end
function prob_detect(prob::SimulatedSystemDetectionProbs{T}, pl_id::Integer)  where T<:SystemDetectionProbsTrait
  if ! (1<=pl_id<=size(prob.pairwise,1) )
    println("#ERROR:  pl_id =", pl_id, " prob.pairwise= ", prob.pairwise)
  end
  @assert 1<=pl_id<=size(prob.pairwise,1)
  return prob.pairwise[pl_id,pl_id]
end
function prob_detect_both_planets(prob::SimulatedSystemDetectionProbs{T}, pl_id::Integer, ql_id::Integer) where T<:SystemDetectionProbsTrait
  @assert 1<=pl_id<=size(prob.pairwise,1)
  @assert 1<=ql_id<=size(prob.pairwise,1)
  prob.pairwise[pl_id,ql_id]
end
prob_detect_n_planets(prob::SimulatedSystemDetectionProbs{T}, n::Integer) where T<:SystemDetectionProbsTrait = 1<=n<=length(prob.n_planets) ? prob.n_planets[n] : 0.0

# Compute sky-averaged transit probabilities from a planetary system with known physical properties, assuming a single host star
#function calc_simulated_system_detection_probs(ps::PlanetarySystemSingleStar, prob_det_if_tr::Vector{Float64}; num_samples::Integer = 1, max_tranets_in_sys::Integer = 10, min_detect_prob_to_be_included::Float64 = 0.0, observer_trait::Type=SkyAveraged)
function calc_simulated_system_detection_probs(ps::PlanetarySystem{StarT}, prob_det_if_tr::Vector{Float64}; num_samples::Integer = 1, max_tranets_in_sys::Integer = 10, min_detect_prob_to_be_included::Float64 = 0.0, observer_trait::Type=SkyAveraged)  where {StarT<:StarAbstract}
  @assert observer_trait <: SystemDetectionProbsTrait
  @assert num_planets(ps) == length(prob_det_if_tr)
  idx_detectable = findall(x->x>0.0,prob_det_if_tr)
  n = length(idx_detectable)
  @assert n <= max_tranets_in_sys       # Make sure memory to store these
  #if n==0
  #   println("# WARNING found no detectable planets based on ",prob_det_if_tr)
  #end
  invalid_prob_flag = false

  #ps_detectable = PlanetarySystemSingleStar(ps,idx_detectable)
  ps_detectable = PlanetarySystem(ps,idx_detectable)

  combo_sample_probs = rand(num_samples)
  combo_cum_probs = zeros(num_samples)
  #sdp = SimulatedSystemDetectionProbs{observer_trait}(prob_det_if_tr[idx_detectable])
  sdp = SimulatedSystemDetectionProbs(observer_trait,prob_det_if_tr[idx_detectable])

  planet_should_transit = zeros(Cint,n)
  for p in 1:num_planets(sdp)
      sdp.pairwise[p,p] = 0.0
  end
  for ntr in 1:min(n,max_tranets_in_sys)  # Loop over number of planets transiting
      sdp.n_planets[ntr] = 0.0
      for combo in combinations(1:n,ntr)  # Loop over specific combinations of detectable planets
        prob_det_this_combo = 1.0
        for p in combo  # Loop over each planet in this combination of detectable planets
            prob_det_this_combo *= prob_det_if_tr[idx_detectable[p]]
        end

        if prob_det_this_combo < min_detect_prob_to_be_included
           continue
        end

	fill!(planet_should_transit,zero(Cint))
	for i in 1:length(combo)
      	    planet_should_transit[combo[i]] = one(Cint)
	end
	local geo_factor::Float64
	if observer_trait == SkyAveraged
	   geo_factor = prob_combo_transits_obs_ave(ps_detectable,planet_should_transit)
	elseif observer_trait == OneObserver
	   geo_factor = prob_combo_transits_one_obs(ps_detectable,planet_should_transit)
	else
	   error(string("typeof(",observer_trait,") is not a valid trait."))
	end
        prob_det_this_combo *= geo_factor

	# Store samples of combinations of planets detected drawn from the full joint multi-planet density
	for i in 1:num_samples
	   if combo_cum_probs[i] < combo_sample_probs[i] <= combo_cum_probs[i]+prob_det_this_combo
	      sdp.combo_detected[i] = combo
	   end
	end
	combo_cum_probs .+= prob_det_this_combo

        sdp.n_planets[ntr] += prob_det_this_combo   # Accumulate the probability of detecting any n planets

        for p in combo                # Accumulate the probability of detecting each planet individually
            sdp.pairwise[p,p] += prob_det_this_combo
        end

        #=
        for pq in combinations(combo,2)                # Accumulate the probability of detecting each planet pair # TODO: OPT: replace with simply calculating integers for pairs to avoid allocations of small arrays
           sdp.pairwise[pq[1],pq[2]] = prob_det_this_combo
           sdp.pairwise[pq[2],pq[1]] = prob_det_this_combo   # TODO OPT: Remove if use symmetric matrix type.
        end
        =#
        if length(combo)>=2
          for pi in 2:length(combo)                # Accumulate the probability of detecting each planet pair # TODO: OPT: replace with simply calculating integers for pairs to avoid allocations of small arrays
            p = combo[pi]
            for qi in 1:(pi-1)
               q = combo[qi]
               sdp.pairwise[p,q] = prob_det_this_combo
               sdp.pairwise[q,p] = prob_det_this_combo   # TODO OPT: Remove if use symmetric matrix type.
            end # qi
          end # pi
        end # if
      end # combo
  end # ntr

  for p in 1:n
      if sdp.pairwise[p,p] > 1.0
          invalid_prob_flag = true
          println(string("Error! Invalid prob for planet ",p,": ", sdp.pairwise[p,p]))
      end
  end

  if invalid_prob_flag
      println("")
      for ntr in 1:min(n,max_tranets_in_sys)
          for combo in combinations(1:n,ntr)
              fill!(planet_should_transit,zero(Cint))
	      for i in 1:length(combo)
      	          planet_should_transit[combo[i]] = one(Cint)
	      end
              if length(combo) == n
                  geo_factor = prob_combo_transits_obs_ave(ps_detectable,planet_should_transit, print_orbit = true)
              else
                  geo_factor = prob_combo_transits_obs_ave(ps_detectable,planet_should_transit)
              end
              println(string("Geo. factor of ",combo," = ",geo_factor))
          end
      end
      println(string("Det. prob. = ", prob_det_if_tr[idx_detectable]))
      println("")
      #quit()
  end
  return sdp
end

if false # WARNING: Complicated and untested
function combine_system_detection_probs(prob::Vector{SimulatedSystemDetectionProbs{T}}, s1::Integer, s2::Integer) where T # WARNING: Complicated and untested
    npl_s1 = min(num_planets(prob[s1]), max_tranets_in_sys)
    npl_s2 = min(num_planets(prob[s2]), max_tranets_in_sys)
    num_planets_across_systems = npl_s1 + npl_s2
    prob_merged = SimulatedSystemDetectionProbs{T}(num_planets_across_systems)
    # Copy probabilities for detecting planets and planet pairs within one system
    prob_merged.detect_planet_if_transits[1:npl_s1] = prob[s1].detect_planet_if_transits
    prob_merged.pairwise[1:npl_s1,1:npl_s1] = prob[s1].pairwise
    offset = npl_s1
    prob_merged.detect_planet_if_transits[offset+1:offset+npl_s2] = prob[s2].detect_planet_if_transits
    prob_merged.pairwise[offset+1:offset+npl_s2,offset+1:offset+npl_s2] = prob[s2].pairwise
    # Calculate probabilities of detecting pairs of planets in different systems, assuming uncorrelated orientations
    for p1 in 1:npl_s1
        for p2 in 1:npl_s2
	    prob_detect_both = prob_detect(prob[s1],p1)*prob_detect(prob[s2],p2)
	    idx1 = p1
	    idx2 = offset+p2
	    prob_merged.pairwise[idx1,idx2] = prob_detect_both
	    prob_merged.pairwise[idx2,idx1] = prob_detect_both
	end
    end
    # Merge probabilities of detecting n_planets, assuming uncorrelated orientations
    p_zero_pl_s1 = 1.0-sum(prob[s1].n_planets)
    p_zero_pl_s2 = 1.0-sum(prob[s2].n_planets)
    for n in 1:min(num_planets_across_systems,max_tranets_in_sys)
        prob_merged.n_planets[n] = p_zero_pl_s1*prob[s2].n_planets[n] + p_zero_pl_s2*prob[s1].n_planets[n]
	for i in 1:n-1
	   prob_merged.n_planets[n] += prob[s1].n_planets[i]*prob[s2].n_planets[n-i] + prob[s1].n_planets[n-i]*prob[s2].n_planets[i]
	end
    end
    # Combine samples of detected planet combinations
    prob_merged.combo_detected = fill(Array{Int64}(undef,0), min(length(prob[s1].combo_detected), length(prob[s2].combo_detected) ) )
    for i in 1:length(prob_merged.combo_detected)
       prob_merged.combo_detected[i] = vcat( prob[s1].combo_detected, prob[s1].combo_detected+offset )
    end
    return prob_merged
end

function select_subset(prob::SimulatedSystemDetectionProbs{T}, idx::Vector{Int64}) # WARNING: Complicated and untested where {T<:SystemDetectionProbsTrait}
    n = length(idx)
    subset = SimulatedSystemDetectionProbs{T}(n)
    subset.detect_planet_if_transits = prob.detect_planet_if_transits[idx]
    subset.pairwise = ones(n,n)
    subset.n_planets = zeros(max(n,length(prob.n_planets)))
    if 1<=n<=length(subset.n_planets)
        subset.n_planets[n] = 1.0
    end
    subset.combo_detected = Array{Int64,1}[collect(1:n)]
    return subset
end


# ASSUMING: Planetary systems for same target are uncorrelated
# Compute sky-averaged transit probabilities from a target with known physical properties
function calc_simulated_system_detection_probs(t::KeplerTarget, sim_param::SimParam ) # WARNING: Complicated and untested
  max_tranets_in_sys = get_int(param,"max_tranets_in_sys",10)
  min_detect_prob_to_be_included = get(param,"max_tranets_in_sys", 0.0)
  s1 = findfirst(num_planets,t.sys)
  if num_planets(t) == num_planets(t.sys[s1])
    # Target has only one system with planets
    prob_det_if_tr = Float64[calc_ave_prob_detect_if_transit(t, s1, p, sim_param) for p in 1:num_planets(t.sys[s1])]
    return calc_simulated_system_detection_probs(t.sys[s1], prob_det_if_tr, max_tranets_in_sys=max_tranets_in_sys, min_detect_prob_to_be_included=min_detect_prob_to_be_included )
  else
    # Target has multiple systems with planets
    sdp = SkyAveragedSystemDetectionProbs[ SkyAveragedSystemDetectionProbs( min(num_planets(t.sys[s]),max_tranets_in_sys) ) for s in 1:length(t.sys) ]
    num_planets_across_systems = 0
    # Calculate detection probabilities for each system separately
    for s in 1:length(t.sys)
        prob_det_if_tr = Float64[calc_ave_prob_detect_if_transit(t, s, p, sim_param) for p in 1:num_planets(t.sys[s])]
    	sdp[s] = calc_simulated_system_detection_probs(t.sys[s], prob_det_if_tr, max_tranets_in_sys=max_tranets_in_sys, min_detect_prob_to_be_included=min_detect_prob_to_be_included )
    	num_planets_across_systems += num_planets(t.sys[s])
    end
    @assert num_planets_across_systems <= max_tranets_in_sys # Make sure memory to store these  # QUERY: DETAIL: Should we relax?
    # Find system ids for first two systems with planets
    #s1, s2 = find_system_detection_probs_with_planets(::Vector{})
    more_than_two_systems_with_planets = false
    s1 = 0
    s2 = 0
    for s in 1:length(t.sys)
	if num_planets(t.sys[s])>=1
	   if s1==0
	      s1 = s
	   elseif s2==0
	      s2 = s
	   else
	      more_than_two_systems_with_planets = true
	   end
	end
    end
    @assert (s1!=0) && (s2!=0)
    @assert !more_than_two_systems_with_planets
    sdp_merged = combine_system_detection_probs(sdp,s1,s2)   # Merge SystemDetectionProbs across systems with common target
    return sdp_merged
  end
end
end


mutable struct ObservedSystemDetectionProbs <: SystemDetectionProbsAbstract          # TODO OPT:  For observed systems (or simulations of observed systems) were we can't know everything.  Is this even used?  Or should we just compute these on the fly, rather than storing them? Do we even want to keep this?
  planet_transits::Vector{Float64}                         # Probability that each planet transits individually for one observer based on actual i, e, and omega
  detect_planet_if_transits::Vector{Float64}               # Probability of detecting each planet given that it transits. Assumes one observer based on actual i, e and omega
  # snr::Vector{Float64}                      # Dimensionless SNR of detection for each planet QUERY: Should we store this here?
end
ObservedSystemDetectionProbs(p::Vector{Float64}) = ObservedSystemDetectionProbs( ones(length(p)), p )
ObservedSystemDetectionProbs(n::Integer) = ObservedSystemDetectionProbs( ones(n), zeros(n) )
ObservedSystemDetectionProbsEmpty() = ObservedSystemDetectionProbs(0)

# Functions common to various types of SystemDetectionProbs
num_planets(prob::ObservedSystemDetectionProbs) = length(prob.detect_planet_if_transits)
prob_detect_if_transits(prob::ObservedSystemDetectionProbs, pl_id::Integer) = prob.detect_planet_if_transits[pl_id]
prob_detect(prob::ObservedSystemDetectionProbs, pl_id::Integer) = prob.planet_transits[pl_id]*prob.detect_planet_if_transits[pl_id]
#prob_detect_both_planets(prob::ObservedSystemDetectionProbs, pl_id::Integer, ql_id::Integer) = prob_detect(prob,pl_id) * prob_detect(prob,ql_id)   # WARNING: Assumes independent.  Intent is for testing CORBITS.  Or should we delete?

#if false
function prob_detect_n_planets(prob::ObservedSystemDetectionProbs, n::Integer)   # WARNING: Assumes independent.  Intent is for testing CORBITS.  May wnat to comment out to prevent accidental use.
  if n<1 || n > num_planets(prob)  return 0.0 end
  sum_prob = 0.0
  for combo in combinations(1:num_planets(prob), n)
      prob_this_combo = 1.0
      for pl_id in combo
          prob_this_combo *= prob_detect(prob,pl_id)
      end
      sum_prob += prob_this_combo
  end
  return sum_prob
end



#=
# Compute transit probabilities for a single observer from a target with known physical properties
function calc_observed_system_detection_probs(targ::KeplerTarget, sim_param::SimParam)
  n = num_planets(targ)
  pdet = zeros(n)
  ptr = zeros(n)
  pl = 1
  for s in 1:length(targ.sys)
      for p in 1:length(targ.sys[s].planet)
         pdet[pl] = calc_prob_detect_if_transit_with_actual_b(targ, s, p, sim_param)
         ptr[pl]  = calc_transit_prob_single_planet_one_obs(targ, s, p)
         pl += 1
      end
  end
  ObservedSystemDetectionProbs( ptr, pdet )
end
=#

# Estimate transit probabilities for a single observer from a target with known physical properties.
if false    # Do we actually want this for anything?
function calc_observed_system_detection_probs(kto::KeplerTargetObs, sim_param::SimParam)
  n = num_planets(kto)
  pdet = ones(n)   # WARNING: We assume all observed objects were detected and we don't have enough info to calcualte a detection probability.  Do we want to do something different?
  ptr = zeros(n)
  for pl in 1:n
      ptr[pl]  = calc_transit_prob_single_planet_approx(kto.obs[pl].period, kto.star.radius, kto.star.mass )
  end
  ObservedSystemDetectionProbs( ptr, pdet )
end
end
