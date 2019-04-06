## ExoplanetsSysSim/src/transit_observations.jl
## (c) 2015 Eric B. Ford

#using Distributions
#include("constants.jl")
include("newPDST.jl")		#includes a function to calculate if given durations match those observed by Kepler for a given period

#  Starting Section of Observables that are actually used
struct TransitPlanetObs
  # ephem::ephemeris_type     # For now hardcode P and t0, see transit_observation_unused.jl to reinstate
  period::Float64             # days
  t0::Float64                 # days
  depth::Float64              # fractional
  duration::Float64           # days; Full-width, half-max-duration until further notice
  # ingress_duration::Float64   # days;  QUERY:  Will we want to use the ingress/egress duration for anything?
end
TransitPlanetObs() = TransitPlanetObs(0.0,0.0,0.0,0.0)

struct StarObs
  radius::Float64      # in Rsol
  mass::Float64        # in Msol
  id::Int64            # row number in stellar dataframe
end


period(obs::TransitPlanetObs) = obs.period
depth(obs::TransitPlanetObs) = obs.depth
duration(obs::TransitPlanetObs) = obs.duration

semimajor_axis(P::Float64, M::Float64) = (grav_const/(4pi^2)*M*P*P)^(1/3)

function semimajor_axis(ps::PlanetarySystemAbstract, id::Integer)
  M = mass(ps.star) + ps.planet[id].mass   # TODO SCI DETAIL: Replace with Jacobi mass?  Not important unless start including TTVs, even then unlikely to matter
  @assert(M>0.0)
  @assert(ps.orbit[id].P>0.0)
  return semimajor_axis(ps.orbit[id].P,M)
end

function calc_transit_depth(t::KeplerTarget, s::Integer, p::Integer)  # WARNING: IMPORTANT: Assumes non-grazing transit
  radius_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  # depth = radius_ratio^2 # TODO SCI DETAIL: Include limb darkening?
  #b = calc_impact_parameter(t.sys[s].planet, p) # If limb darkening should know about which chord the planet takes set b to impact parameter, rather than 0.0.
  #depth = depth_small_planet(radius_ratio, b, t.sys[s].star.ld)   # TODO Add limb darkening here .  When limb darkening is implemented, then we should update calc_snr_if_transit in transit_detection_model to not convert transit depth to tps deth.
  depth = depth_at_midpoint(radius_ratio, t.sys[s].star.ld)   # TODO Add limb darkening here.  When limb darkening is implemented, then we should update calc_snr_if_transit in transit_detection_model to not convert transit depth to tps deth.
  depth *=  flux(t.sys[s].star)/flux(t)                      # Flux ratio accounts for dilution
end

function calc_transit_duration_central_circ_small_angle_approx(ps::PlanetarySystemAbstract, pl::Integer)
  duration = rsol_in_au*ps.star.radius * ps.orbit[pl].P /(pi*semimajor_axis(ps,pl) )
end

calc_transit_duration_central_circ_small_angle_approx(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_circ_small_angle_approx(t.sys[s],p)

function calc_transit_duration_central_circ_with_arcsin(ps::PlanetarySystemAbstract, pl::Integer)
  asin_arg = rsol_in_au*ps.star.radius/semimajor_axis(ps,pl)
  duration = ps.orbit[pl].P/pi * (asin_arg < 1.0 ? asin(asin_arg) : 1.0)
end
calc_transit_duration_central_circ_with_arcsin(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_circ_with_arcsin(t.sys[s],p)


#calc_transit_duration_central_circ(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_central_circ_small_angle_approx(ps,pl)
calc_transit_duration_central_circ(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_central_circ_with_arcsin(ps,pl)

calc_transit_duration_central_circ(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_circ(t.sys[s],p)


function calc_transit_duration_central_small_angle_approx(ps::PlanetarySystemAbstract, pl::Integer)
  ecc = ps.orbit[pl].ecc
  sqrt_one_minus_ecc_sq = sqrt((1+ecc)*(1-ecc))
  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  vel_fac = sqrt_one_minus_ecc_sq/one_plus_e_sin_w
  duration = calc_transit_duration_central_circ_small_angle_approx(ps,pl) * vel_fac
end
calc_transit_duration_central_small_angle_approx(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_small_angle_approx(t.sys[s],p)

# It seems above should be good enough, but we can try one of the following just to eliminate potential approximation errors.
function calc_transit_duration_central_winn2010(ps::PlanetarySystemAbstract, pl::Integer)
  ecc = ps.orbit[pl].ecc
  sqrt_one_minus_ecc_sq = sqrt((1+ecc)*(1-ecc))
  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  vel_fac = sqrt_one_minus_ecc_sq/one_plus_e_sin_w
  radial_separation_over_a = (1+ecc)*(1-ecc)/one_plus_e_sin_w
  asin_arg = rsol_in_au*ps.star.radius/(semimajor_axis(ps,pl))
  # Based on Winn 2010
  duration = ( asin_arg<1.0 ?  asin(asin_arg) : 1.0 ) * ps.orbit[pl].P*radial_separation_over_a/(pi*sqrt_one_minus_ecc_sq)
end
calc_transit_duration_central_winn2010(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_winn2010(t.sys[s],p)

function calc_transit_duration_central_kipping2010(ps::PlanetarySystemAbstract, pl::Integer)
  ecc = ps.orbit[pl].ecc
  sqrt_one_minus_ecc_sq = sqrt((1+ecc)*(1-ecc))
  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  vel_fac = sqrt_one_minus_ecc_sq/one_plus_e_sin_w
  radial_separation_over_a = (1+ecc)*(1-ecc)/one_plus_e_sin_w
  # Based on pasting cos i = 0 into Eqn 15 from Kipping 2010
  asin_arg = rsol_in_au*ps.star.radius/(semimajor_axis(ps,pl)* radial_separation_over_a)
  duration = ps.orbit[pl].P*radial_separation_over_a^2/(pi*sqrt_one_minus_ecc_sq) * ( asin_arg<1.0 ?  asin(asin_arg) : 1.0 )
end
calc_transit_duration_central_kipping2010(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central_kipping2010(t.sys[s],p)

#calc_transit_duration_central(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_central_small_angle_approx(ps,pl)
#calc_transit_duration_central(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_central_winn2010(ps,pl)
calc_transit_duration_central(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_central_kipping2010(ps,pl)

calc_transit_duration_central(t::KeplerTarget, s::Integer, p::Integer) = calc_transit_duration_central(t.sys[s],p)

function calc_transit_duration_factor_for_impact_parameter_b(b::T, p::T)  where T <:Real
  @assert(zero(b)<=b)         # b = Impact Parameter
  @assert(zero(p)<=p<one(p))  # p = R_p/R_star
  if b < 1-3p                 # Far enough from grazing for approximation
        #duration_ratio = sqrt(1-b^2)  # Approximation to (sqrt((1+p)^2-b^2)+sqrt((1-p)^2-b^2))/2, which is itself an approximation
        duration_ratio = sqrt((1-b)*(1+b))  # Approximation to (sqrt((1+p)^2-b^2)+sqrt((1-p)^2-b^2))/2, which is itself an approximation
    elseif b < 1-p            # Planet is fully inscribed at mid-transit
        #duration_ratio = (sqrt((1+p)^2-b^2)+sqrt((1-p)^2-b^2))/2 # Average of full and flat transit durations approximates to duration between center of planet being over limb of star
        duration_ratio = (sqrt(((1+p)+b)*((1+p)-b))+sqrt(((1-p)+b)*((1-p)-b)))/2 # Average of full and flat transit durations approximates to duration between center of planet being over limb of star
    elseif b < 1+p            # Planet never fully inscribed by star
        #duration_ratio = sqrt((1+p)^2-b^2)/2  # /2 since now triangular
        duration_ratio = sqrt(((1+p)+b)*((1+p)-b))/2  # /2 since now triangular
    else                      # There's no transit
        duration_ratio = zero(b)
    end
    return duration_ratio
end


function calc_effective_transit_duration_factor_for_impact_parameter_b(b::T, p::T)  where T <:Real
  @assert(zero(b)<=b)         # b = Impact Parameter
  @assert(zero(p)<=p<one(p))  # p = R_p/R_star
  if b < 1-3p                 # Far enough from grazing for approximation
        duration_ratio = sqrt((1+b)*(1-b))  # Approximation to (sqrt((1+p)^2-b^2)+sqrt((1-p)^2-b^2))/2, which is itself an approximation
        area_ratio = one(p)
    elseif b < 1-p            # Planet is fully inscribed at mid-transit
        #duration_ratio = (sqrt((1+p)^2-b^2)+sqrt((1-p)^2-b^2))/2 # Average of full and flat transit durations approximates to duration between center of planet being over limb of star
        duration_ratio = (sqrt(((1+p)+b)*((1+p)-b))+sqrt(((1-p)+b)*((1-p)-b)))/2 # Average of full and flat transit durations approximates to duration between center of planet being over limb of star
        area_ratio = one(p)
    elseif b < 1+p            # Planet never fully inscribed by star
        #duration_ratio = sqrt((1+p)^2-b^2)/2  # /2 since now triangular
        duration_ratio = sqrt(((1+p)+b)*((1+p)-b))/2  # /2 since now triangular
        #area_ratio = (p^2*acos((b^2+p^2-1)/(2*b*p))+acos((b^2+1-p^2)/(2b))-0.5*sqrt((1+p-b)*(p+b-1)*(1-p+b)*(1+p+b))) / (pi*p^2)
        #area_ratio = (p^2*acos((b^2-(1-p)*(1+p))/(2*b*p))+acos((b^2+(1+p)*(1-p))/(2b))-0.5*sqrt((1+p-b)*(p+b-1)*(1-p+b)*(1+p+b))) / (pi*p^2)
        acos_arg1 = max(-1.0,min(1.0,(b^2-(1-p)*(1+p))/(2*b*p)))
        acos_arg2 = max(-1.0,min(1.0,(b^2+(1+p)*(1-p))/(2b)))
        sqrt_arg = max(0.0,(1+p-b)*(p+b-1)*(1-p+b)*(1+p+b))
        area_ratio = (p^2*acos(acos_arg1)+acos(acos_arg2)-0.5*sqrt(sqrt_arg)) / (pi*p^2)
    else                      # There's no transit
        duration_ratio = zero(b)
        area_ratio = zero(p)
    end
    return duration_ratio*area_ratio
end


# How SNR is affected for grazing transits due to not all of planet blocking starlight at mid-transit.
# Assumes uniform surface brightness star
# Expression comes from Eqn 14 of http://mathworld.wolfram.com/Circle-CircleIntersection.html
function calc_depth_correction_for_grazing_transit(b::T, p::T)  where T <:Real
  @assert(zero(b)<=b)         # b = Impact Parameter
  @assert(zero(p)<=p<one(p))  # p = R_p/R_star
  if b < 1-p                  # Planet fully inscribed by star
        area_ratio = one(p)
    elseif b < 1+p            # Planet never fully inscribed by star
        #area_ratio = (p^2*acos((b^2+p^2-1)/(2*b*p))+acos((b^2+1-p^2)/(2b))-0.5*sqrt((1+p-b)*(p+b-1)*(1-p+b)*(1+p+b))) / (pi*p^2)
        acos_arg1 = max(-1.0,min(1.0,(b^2-(1-p)*(1+p))/(2*b*p)))
        acos_arg2 = max(-1.0,min(1.0,(b^2+(1+p)*(1-p))/(2b)))
        sqrt_arg = max(0.0,(1+p-b)*(p+b-1)*(1-p+b)*(1+p+b))
        area_ratio = (p^2*acos(acos_arg1)+acos(acos_arg2)-0.5*sqrt(sqrt_arg)) / (pi*p^2)
    else                      # There's no transit
        area_ratio = zero(p)
    end
    return area_ratio
end


function calc_transit_duration_small_angle_approx(ps::PlanetarySystemAbstract, pl::Integer)
  a = semimajor_axis(ps,pl)
  @assert a>=zero(a)
  ecc = ps.orbit[pl].ecc
  @assert zero(ecc)<=ecc<=one(ecc)
  b = calc_impact_parameter(ps, pl)
  size_ratio = ps.planet[pl].radius/ps.star.radius
  @assert !isnan(b)
  @assert zero(b)<=b
  if b>one(b)+size_ratio
     return zero(b)
  end
  duration_central_circ = calc_transit_duration_central_circ(ps,pl)
  duration_ratio_for_impact_parameter = calc_transit_duration_factor_for_impact_parameter_b(b,size_ratio)

  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  sqrt_one_minus_e_sq = sqrt((1+ecc)*(1-ecc))
  vel_fac = sqrt_one_minus_e_sq / one_plus_e_sin_w

  duration = duration_central_circ * duration_ratio_for_impact_parameter * vel_fac
end
calc_transit_duration_small_angle_approx(t::KeplerTarget, s::Integer, p::Integer ) = calc_transit_duration_small_angle_approx(t.sys[s],p)

function calc_transit_duration_winn2010(ps::PlanetarySystemAbstract, pl::Integer)
  a = semimajor_axis(ps,pl)
  @assert a>=zero(a)
  ecc = ps.orbit[pl].ecc
  @assert zero(ecc)<=ecc<=one(ecc)
  b = calc_impact_parameter(ps, pl)
  size_ratio = ps.planet[pl].radius/ps.star.radius
  @assert !isnan(b)
  @assert zero(b)<=b
  if b>one(b)+size_ratio
     return zero(b)
  end
  duration_central_circ = calc_transit_duration_central_circ(ps,pl)
  arcsin_circ_central = pi/ps.orbit[pl].P*duration_central_circ

  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  sqrt_one_minus_e_sq = sqrt((1+ecc)*(1-ecc))
  vel_fac = sqrt_one_minus_e_sq / one_plus_e_sin_w
  radial_separation_over_a = (1+ecc)*(1-ecc)/one_plus_e_sin_w
  duration_ratio_for_impact_parameter = calc_transit_duration_factor_for_impact_parameter_b(b,size_ratio)

  # WARNING: This is technically an approximation.  It avoids small angle for non-grazing transits, but does use a variant of the small angle approximation for nearly and grazing transits.
  asin_arg = (arcsin_circ_central * duration_ratio_for_impact_parameter)
  duration = ps.orbit[pl].P/pi * radial_separation_over_a/sqrt_one_minus_e_sq * (asin_arg < 1.0 ? asin(asin_arg) : 1.0)
  duration = duration_central_cric * radial_separation_over_a/sqrt_one_minus_e_sq * (asin_arg < 1.0 ? asin(asin_arg) : 1.0)/arcsin_circ_central
end
calc_transit_duration_winn2010(t::KeplerTarget, s::Integer, p::Integer ) = calc_transit_duration_winn2010(t.sys[s],p)

function calc_transit_duration_kipping2010(ps::PlanetarySystemAbstract, pl::Integer)
  a = semimajor_axis(ps,pl)
  @assert a>=zero(a)
  ecc = ps.orbit[pl].ecc
  @assert zero(ecc)<=ecc<=one(ecc)
  b = calc_impact_parameter(ps, pl)
  size_ratio = ps.planet[pl].radius/ps.star.radius
  @assert !isnan(b)
  @assert zero(b)<=b
  if b>one(b)+size_ratio
     return zero(b)
  end
  duration_central_circ = calc_transit_duration_central_circ(ps,pl)
  arcsin_circ_central = pi/ps.orbit[pl].P*duration_central_circ

  one_plus_e_sin_w = 1+ecc*sin(ps.orbit[pl].omega)
  sqrt_one_minus_e_sq = sqrt((1+ecc)*(1-ecc))
  vel_fac = sqrt_one_minus_e_sq / one_plus_e_sin_w
  radial_separation_over_a = (1+ecc)*(1-ecc)/one_plus_e_sin_w
  duration_ratio_for_impact_parameter = calc_transit_duration_factor_for_impact_parameter_b(b,size_ratio)

  # WARNING: This is technically an approximation (see Kipping 2010 Eqn 15).  It avoids small angle for non-grazing transits, but does use a variant of the small angle approximation for nearly and grazing transits.
  asin_arg = (arcsin_circ_central * duration_ratio_for_impact_parameter/radial_separation_over_a)
  duration = ps.orbit[pl].P/pi * radial_separation_over_a^2/sqrt_one_minus_e_sq * (asin_arg < 1.0 ? asin(asin_arg) : 1.0)
  #duration = duration_central_circ * radial_separation_over_a^2/sqrt_one_minus_e_sq * (asin_arg < 1.0 ? asin(asin_arg) : 1.0)/arcsin_circ_central
end
calc_transit_duration_kipping2010(t::KeplerTarget, s::Integer, p::Integer ) = calc_transit_duration_kipping2010(t.sys[s],p)

#calc_transit_duration(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_small_angle_approx(ps,pl)
#calc_transit_duration(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_winn2010(ps,pl)
calc_transit_duration(ps::PlanetarySystemAbstract, pl::Integer) = calc_transit_duration_kipping2010(ps,pl)
calc_transit_duration(t::KeplerTarget, s::Integer, p::Integer ) = calc_transit_duration(t.sys[s],p)

function calc_expected_num_transits(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
 period = t.sys[s].orbit[p].P
 exp_num_transits = t.duty_cycle * t.data_span/period
 return exp_num_transits
end

include("transit_detection_model.jl")
include("transit_prob_geometric.jl")

const has_sc_bit_array_size = 7*8        # WARNING: Must be big enough given value of num_quarters (assumed to be <=17)

mutable struct KeplerTargetObs                        # QUERY:  Do we want to make this type depend on whether the catalog is based on simulated or real data?
  obs::Vector{TransitPlanetObs}
  sigma::Vector{TransitPlanetObs}           # Simplistic approach to uncertainties for now.  QUERY: Should estimated uncertainties be part of Observations type?
  # phys_id::Vector{Tuple{Int32,Int32}}     # So we can lookup the system's properties # Commented out since Not used

  prob_detect::SystemDetectionProbsAbstract  # QUERY: Specialize type of prob_detect depending on whether for simulated or real data?

  #has_sc::Vector{Bool}                      # TODO OPT: Make Immutable Vector or BitArray to reduce memory use?
  has_sc::BitArray{1}                        # WARNING: Changed from Array{Bool}.  Alternatively, we could try StaticArray{Bool} so fixed size?  Do we even need to keep this?

  star::StarObs
end
#KeplerTargetObs(n::Integer) = KeplerTargetObs( fill(TransitPlanetObs(),n), fill(TransitPlanetObs(),n), fill(tuple(0,0),n),  ObservedSystemDetectionProbsEmpty(),  fill(false,num_quarters), StarObs(0.0,0.0) )
KeplerTargetObs(n::Integer) = KeplerTargetObs( fill(TransitPlanetObs(),n), fill(TransitPlanetObs(),n), ObservedSystemDetectionProbsEmpty(),  falses(has_sc_bit_array_size), StarObs(0.0,0.0,0) )
num_planets(t::KeplerTargetObs) = length(t.obs)

function calc_target_obs_sky_ave(t::KeplerTarget, sim_param::SimParam)
   max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
   transit_noise_model = get_function(sim_param,"transit_noise_model")
   min_detect_prob_to_be_included = 0.0  # get_real(sim_param,"min_detect_prob_to_be_included")
   num_observer_samples = 1 # get_int(sim_param,"num_viewing_geometry_samples")
  vetting_efficiency = get_function(sim_param,"vetting_efficiency")

  np = num_planets(t)
  obs = Array{TransitPlanetObs}(undef,np)
  sigma = Array{TransitPlanetObs}(undef,np)
  #id = Array{Tuple{Int32,Int32}}(np)
  #id = Array{Tuple{Int32,Int32}}(np)
  ns = length(t.sys)
  sdp_sys = Array{SystemDetectionProbsAbstract}(undef,ns)
  i = 1
  for (s,sys) in enumerate(t.sys)
    pdet = zeros(num_planets(sys))
    for (p,planet) in enumerate(sys.planet)
      if get(sim_param,"verbose",false)
         println("# s=",s, " p=",p," num_sys= ",length(t.sys), " num_pl= ",num_planets(sys) )
      end
        period = sys.orbit[p].P
        duration_central = calc_transit_duration_central(t,s,p)
        size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
        depth = calc_transit_depth(t,s,p)
        ntr = calc_expected_num_transits(t, s, p, sim_param)

        # cdpp_central = interpolate_cdpp_to_duration(t, duration_central)
        # snr_central = calc_snr_if_transit_cdpp(t, depth, duration_central, cdpp_central, sim_param, num_transit=ntr)
        # pdet_ave = calc_ave_prob_detect_if_transit_from_snr_cdpp(t, snr_central, period, duration_central, size_ratio, cdpp_central, sim_param, num_transit=ntr)

        kepid = StellarTable.star_table(t.sys[s].star.id, :kepid)
        osd_duration_central = get_legal_durations(period,duration_central)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
        osd_central = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration_central)
        if osd_duration_central > duration_central				#use a correcting factor if this duration is lower than the minimum searched for this period.
	   osd_central = osd_central*osd_duration_central/duration_central
	end
        snr_central = calc_snr_if_transit(t, depth, duration_central, osd_central, sim_param, num_transit=ntr)
        pdet_ave = calc_ave_prob_detect_if_transit_from_snr(t, snr_central, period, duration_central, size_ratio, osd_central, sim_param, num_transit=ntr)

	add_to_catalog = pdet_ave > min_detect_prob_to_be_included  # Include all planets with sufficient detection probability

	if add_to_catalog
           pdet_central = calc_prob_detect_if_transit(t, snr_central, period, duration_central, sim_param, num_transit=ntr)
           threshold_pdet_ratio = rand()
	   hard_max_num_b_tries = 100
	   max_num_b_tries = min_detect_prob_to_be_included == 0. ? hard_max_num_b_tries : min(hard_max_num_b_tries,convert(Int64,1/min_detect_prob_to_be_included))
           # We compute measurement noise based on a single value of b.  We draw from a uniform distribution for b and then using rejection sampling to reduce probability of higher impact parameters
           pdet_this_b = 0.0
           for j in 1:max_num_b_tries
              b = rand()
              transit_duration_factor = calc_effective_transit_duration_factor_for_impact_parameter_b(b,size_ratio)

	      duration = duration_central * transit_duration_factor   # WARNING:  Technically, this duration may be slightly reduced for grazing cases to account for reduction in SNR due to planet not being completely inscribed by star at mid-transit.  But this will be a smaller effect than limb-darkening for grazing transits.  Also, makes a variant of the small angle approximation

              # cdpp = interpolate_cdpp_to_duration(t, duration)
              # snr = snr_central * (cdpp_central/cdpp) * sqrt(transit_duration_factor)

              osd_duration = get_legal_durations(period,duration)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
              osd = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration)
              if osd_duration > duration				#use a correcting factor if this duration is lower than the minimum searched for this period.
	          osd = osd*osd_duration/duration
	      end
              snr = snr_central * (osd_central/osd)

              pdet_this_b = calc_prob_detect_if_transit(t, snr, period, duration, sim_param, num_transit=ntr)
              pvet = vetting_efficiency(t.sys[s].planet[p].radius, period) # TODO: Ask Danely about this line

              if pdet_this_b >= threshold_pdet_ratio * pdet_central
                  #println("# Adding pdet_this_b = ", pdet_this_b, " pdet_c = ", pdet_central, " snr= ",snr, " cdpp= ",cdpp, " duration= ",duration, " b=",b, " u01= ", threshold_pdet_ratio)
	         pdet[p] = pdet_ave*pvet
                 obs[i], sigma[i] = transit_noise_model(t, s, p, depth, duration, snr, ntr, b=b)
      	         i += 1
                 break
              end # if add to obs and sigma lists
           end # for j
	else # add_to_catalog
	   # Do anything for planets that are extremely unlikely to be detected even if they were to transit?
	end
    end
    resize!(obs,i-1)
    resize!(sigma,i-1)
    sdp_sys[s] = calc_simulated_system_detection_probs(sys, pdet, max_tranets_in_sys=max_tranets_in_sys, min_detect_prob_to_be_included=min_detect_prob_to_be_included, num_samples=num_observer_samples)
  end
  # TODO SCI DETAIL: Combine sdp_sys to allow for target to have multiple planetary systems
  s1 = findfirst(x->num_planets(x)>0,sdp_sys)  # WARNING IMPORTANT: For now just take first system with planets
  if s1 == nothing
     s1 = 1
  end
  sdp_target = sdp_sys[s1]

  has_no_sc = falses(has_sc_bit_array_size)
  star_obs = StarObs( t.sys[1].star.radius, t.sys[1].star.mass, t.sys[1].star.id )  # NOTE:  This sets the observed star properties to be those in the stellar catalog.  If want to incorporate uncertainty in stellar properties, that would be done elsewhere when translating depths into planet radii.
  #return KeplerTargetObs(obs, sigma, id, sdp_target, has_no_sc, star_obs )
  return KeplerTargetObs(obs, sigma, sdp_target, has_no_sc, star_obs )
end


function calc_target_obs_single_obs(t::KeplerTarget, sim_param::SimParam)
  # max_tranets_in_sys = get_int(sim_param,"max_tranets_in_sys")
   transit_noise_model = get_function(sim_param,"transit_noise_model")
   min_detect_prob_to_be_included = 0.0  # get_real(sim_param,"min_detect_prob_to_be_included")
  transit_noise_model = get_function(sim_param,"transit_noise_model")

  np = num_planets(t)
  obs = Array{TransitPlanetObs}(undef,np)
  sigma = Array{TransitPlanetObs}(undef,np)
  ns = length(t.sys)
  sdp_sys = Array{ObservedSystemDetectionProbs}(undef,ns)
  i = 1
  cuantos = 1000			#indicator for testing OSD interpolator.
  for (s,sys) in enumerate(t.sys)
    pdet = zeros(num_planets(sys))
    for (p,planet) in enumerate(sys.planet)
      if get(sim_param,"verbose",false)
         println("# s=",s, " p=",p," num_sys= ",length(t.sys), " num_pl= ",num_planets(sys) )
      end
        duration = calc_transit_duration(t,s,p)
	if duration <= 0.
	   continue
	end
        period = sys.orbit[p].P
        ntr = calc_expected_num_transits(t, s, p, sim_param)
        depth = calc_transit_depth(t,s,p)
        # Apply correction to snr if grazing transit
        size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
        b = calc_impact_parameter(t.sys[s],p)
        snr_correction = calc_depth_correction_for_grazing_transit(b,size_ratio)
        depth *= snr_correction

        # cdpp = interpolate_cdpp_to_duration(t, duration)
        # snr = calc_snr_if_transit_cdpp(t, depth, duration, cdpp, sim_param, num_transit=ntr)

        kepid = StellarTable.star_table(t.sys[s].star.id, :kepid)
        osd_duration = get_legal_durations(period,duration)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
        osd = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration)
        if osd_duration > duration				#use a correcting factor if this duration is lower than the minimum searched for this period.
	   osd = osd*osd_duration/duration
        end
        snr = calc_snr_if_transit(t, depth, duration, osd, sim_param, num_transit=ntr)

        pdet[p] = calc_prob_detect_if_transit(t, snr, period, duration, sim_param, num_transit=ntr)

	if pdet[p] > min_detect_prob_to_be_included
           pvet = vetting_efficiency(t.sys[s].planet[p].radius, period) # TODO: Ask Danely about this line
           pdet[p] *= pvet
            obs[i], sigma[i] = transit_noise_model(t, s, p, depth, duration, snr, ntr)
      	   i += 1
	end
    end
    resize!(obs,i-1)
    resize!(sigma,i-1)
    sdp_sys[s] = ObservedSystemDetectionProbs(pdet)
  end
  # TODO SCI DETAIL: Combine sdp_sys to allow for target to have multiple planetary systems
  s1 = findfirst(x->num_planets(x)>0,sdp_sys)  # WARNING: For now just take first system with planets, assumes not two stars wht planets in one target
  if s1 == nothing
     s1 = 1
  end
  sdp_target = sdp_sys[s1]

  has_no_sc = falses(3*num_quarters)
  star_obs = StarObs( t.sys[1].star.radius, t.sys[1].star.mass, t.sys[1].star.id )  # NOTE: This just copies star properties directly
  return KeplerTargetObs(obs, sigma, sdp_target, has_no_sc, star_obs )
end


function test_transit_observations(sim_param::SimParam; verbose::Bool=false)  # TODO TEST: Add more tests
  #transit_param = TransitParameter( EphemerisLinear(10.0, 0.0), TransitShape(0.01, 3.0/24.0, 0.5) )
  generate_kepler_target = get_function(sim_param,"generate_kepler_target")
   max_it = 100000
  local obs
  for i in 1:max_it
    target = generate_kepler_target(sim_param)::KeplerTarget
    while num_planets(target) == 0
      target = generate_kepler_target(sim_param)::KeplerTarget
    end
    #calc_transit_prob_single_planet_one_obs(target,1,1)
    calc_transit_prob_single_planet_obs_ave(target,1,1)
    obs = calc_target_obs_single_obs(target,sim_param)
    obs = calc_target_obs_sky_ave(target,sim_param)
    if verbose && (num_planets(obs) > 0)
      println("# i= ",string(i)," np= ",num_planets(obs), " obs= ", obs )
      break
    end
  end
  return obs
end


randtn() = rand(TruncatedNormal(0.0,1.0,-0.999,0.999))

function transit_noise_model_no_noise(t::KeplerTarget, s::Integer, p::Integer, depth::Float64, duration::Float64, snr::Float64, num_tr::Float64; b::Float64 = 0.0)
  period = t.sys[s].orbit[p].P
  #t0 = rand(Uniform(0.0,period))    # WARNING: Not being calculated from orbit
  t0 = period*rand()    # WARNING: Not being calculated from orbit
  sigma_period = 0.0
  sigma_t0 = 0.0
  sigma_depth =  0.0
  sigma_duration =  0.0
  sigma = TransitPlanetObs( sigma_period, sigma_t0, sigma_depth, sigma_duration )
  obs = TransitPlanetObs( period, t0, depth,duration)
  return obs, sigma
end

function transit_noise_model_fixed_noise(t::KeplerTarget, s::Integer, p::Integer, depth::Float64, duration::Float64, snr::Float64, num_tr::Float64; b::Float64 = 0.0)
  period = t.sys[s].orbit[p].P
  #t0 = rand(Uniform(0.0,period))   # WARNING: Not being calculated from orbit
  t0 = period*rand()    # WARNING: Not being calculated from orbit

  sigma_period = 1e-6
  sigma_t0 = 1e-4
  sigma_depth =  0.1
  sigma_duration =  0.01

  sigma = TransitPlanetObs( sigma_period, sigma_t0, sigma_depth, sigma_duration)
  #obs = TransitPlanetObs( period, t0, depth, duration)
  obs = TransitPlanetObs( period*(1.0+sigma.period*randtn()), t0*(1.0+sigma.period*randtn()), depth*(1.0+sigma.depth*randtn()),duration*(1.0+sigma.duration*randtn()))
  return obs, sigma
end

function transit_noise_model_diagonal(t::KeplerTarget, s::Integer, p::Integer, depth::Float64, duration::Float64, snr::Float64, num_tr::Float64; b::Float64 = calc_impact_parameter(t.sys[s],p) )
  period = t.sys[s].orbit[p].P
  #t0 = rand(Uniform(0.0,period))    # WARNING: Not being calculated from orbit
  t0 = period*rand()    # WARNING: Not being calculated from orbit

	# Use variable names from Price & Rogers
	one_minus_e2 = (1-t.sys[s].orbit[p].ecc)*(1+t.sys[s].orbit[p].ecc)
	a_semimajor_axis = semimajor_axis(t.sys[s],p)

	tau0 = rsol_in_au*t.sys[s].star.radius*period/(a_semimajor_axis*2pi)
	tau0 *= sqrt(one_minus_e2)/(1+t.sys[s].orbit[p].ecc*sin(t.sys[s].orbit[p].omega))
	r = t.sys[s].planet[p].radius/t.sys[s].star.radius
	sqrt_one_minus_b2 = (0.0<=b<1.0) ? sqrt((1-b)*(1+b)) : 0.0
	@assert(sqrt_one_minus_b2>0.0)

 	T = 2*tau0*sqrt_one_minus_b2
	tau = 2*tau0*r/sqrt_one_minus_b2
	Ttot = period
	I = LC_integration_time      # WARNING: Assumes LC only
	Lambda_eff = LC_rate * num_tr # calc_expected_num_transits(t, s, p, sim_param)
	delta = depth
	#sigma = t.cdpp[1,1]
	sigma = interpolate_cdpp_to_duration(t, duration)

	# Price & Rogers Eqn A8 & Table 1 # WARNING: Someone should check Eqns
	tau3 = tau^3
	I3 = I^3
	a1 = (10*tau3+2*I^3-5*tau*I^2)/tau3
	a2 = (5*tau3+I3-5*tau*tau*I)/tau3
	a3 = (9*I^5*Ttot-40*tau3*I*I*Ttot+120*tau^4*I*(3*Ttot-2*tau))/tau^6
	a4 = (a3*tau^5+I^4*(54*tau-35*Ttot)-12*tau*I3*(4*tau+Ttot)+360*tau^4*(tau-Ttot))/tau^5
	a5 = (a2*(24T*T*(I-3*tau)-24*T*Ttot*(I-3*tau))+tau3*a4)/tau3
	a6 = (3*tau*tau+T*(I-3*tau))/(tau*tau)
	a7 = (-60*tau^4+12*a2*tau3*T-9*I^4+8*tau*I3+40*tau3*I)/(tau^4)
	a8 = (2T-Ttot)/tau
	a9 = (-3*tau*tau*I*(-10*T*T+10*T*Ttot+I*(2*I+5*Ttot))-I^4*Ttot+8*tau*I3*Ttot)/(tau^5)
	a10 = ((a9+60)*tau*tau+10*(-9*T*T+9*T*Ttot+I*(3*I+Ttot))-75*tau*Ttot)/(tau*tau)
	a11 = (I*Ttot-3*tau*(Ttot-2*tau))/(tau*tau)
	a12 = (-360*tau^5-24*a2*tau3*T*(I-3*tau)+9*I^5-35*tau*I^4-12*tau*tau*I3-40*tau3*I*I+360*tau^4*I)/(tau^5)
	a13 = (-3*I3*(8*T*T-8*T*Ttot+3*I*Ttot)+120*tau*tau*T*I*(T-Ttot)+8*tau*I3*Ttot)/tau^5
	a14 = (a13*tau*tau+40*(-3*T*T+3*T*Ttot+I*Ttot)-60*tau*Ttot)/(tau*tau)
	a15 = (2*I-6*tau)/tau
	b1  = (6*I*I-3*I*Ttot+tau*Ttot)/(I*I)
	b2  = (tau*T+3*I*(I-T))/(I*I)
	b3 = (tau3-12*T*I*I+8*I3+20*tau*I*I-8*tau*tau*I)/I3
	b4 = (6*T*T-6*T*Ttot+I*(5*Ttot-4*I))/(I*I)
	b5 = (10*I-3*tau)/I
	b6 = (12*b4*I3+4*tau*(-6*T*T+6T*Ttot+I*(13*Ttot-30*I)))/I3
	b7 = (b6*I^5+4*tau*tau*I*I*(12*I-11*Ttot)+tau3*I*(11*Ttot-6*I)-tau^4*Ttot)/I^5
	b8 = (3T*T-3*T*Ttot+I*Ttot)/(I*I)
	b9 = (8*b8*I^4+20*tau*I*I*Ttot-8*tau*tau*I*Ttot+tau3*Ttot)/I^4
	b10 =  (-tau^4+24*T*I*I*(tau-3I)+60*I^4+52*tau*I3-44*tau*tau*I*I+11*tau3*I)/I^4
	b11 =  (-15*b4*I3+10*b8*tau*I*I+15*tau*tau*(2*I-Ttot))/I3
	b12 =  (b11*I^5+2*tau3*I*(4*Ttot-3*I)-tau^4*Ttot)/I^5
	b13 =  (Ttot-2*T)/I
	b14 =  (6*I-2*tau)/I

        Q = snr/sqrt(num_tr)
        sigma_t0 = tau>=I ?  sqrt(0.5*tau*T/(1-I/(3*tau)))/Q : sqrt(0.5*I*T/(1-tau/(3*I)))/Q
        sigma_period = sigma_t0/sqrt(num_tr)
        sigma_duration = tau>=I ? sigma*sqrt(abs(6*tau*a14/(delta*delta*a5)) /Lambda_eff )  : sigma*sqrt(abs(6*I*b9/(delta*delta*b7)) / Lambda_eff)
        sigma_depth = tau>=I ? sigma*sqrt(abs(-24*a11*a2/(tau*a5)) / Lambda_eff)  : sigma*sqrt(abs(24*b1/(I*b7)) / Lambda_eff)

        sigma_obs = TransitPlanetObs( sigma_period, sigma_t0, sigma_depth, sigma_duration )

        local obs
        if true     # Assume uncertainties uncorrelated (Diagonal)
  	    obs = TransitPlanetObs( period*(1.0+sigma_obs.period*randtn()), t0*(1.0+sigma_obs.period*randtn()), depth*(1.0+sigma_obs.depth*randtn()),duration*(1.0+sigma_obs.duration*randtn()))
        else        # TODO SCI DETAIL:  Account for correlated uncertaintties in transit parameters
            cov = zeros(4,4)
        if tau>=I
	# cov[0,0] = -3*tau/(delta*delta*a15)
	cov[1,1] = 24*tau*a10/(delta*delta*a5)
	cov[1,2] = cov[2,1] = 36*a8*tau*a1/(delta*delta*a5)
	cov[1,3] = cov[3,1] = -12*a11*a1/(delta*a5)
	cov[1,4] = cov[4,1] = -12*a6*a1/(delta*a5)
	cov[2,2] = 6*tau*a14/(delta*delta*a5)
	cov[2,3] = cov[3,2] = 72*a8*a2/(delta*a5)
	cov[2,4] = cov[4,2] = 6*a7/(delta*a5)
	cov[3,3] = -24*a11*a2/(tau*a5)
	cov[3,4] = cov[4,3] = -24*a6*a2/(tau*a5)
	cov[4,4] = a12/(tau*a5)
        else
	# cov[0,0] = 3*I/(delta*delta*b14)
	cov[1,1] = -24*I*I*b12/(delta*delta*tau*b7)
	cov[1,2] = cov[2,1] = 36*I*b13*b5/(delta*delta*b7)
	cov[1,3] = cov[3,1] = 12*b5*b1/(delta*b7)
	cov[1,4] = cov[4,1] = 12*b5*b2/(delta*b7)
	cov[2,2] = 6*I*b9/(delta*delta*b7)
	cov[2,3] = cov[3,2] = 72*b13/(delta*b7)
	cov[2,4] = cov[4,2] = 6*b3/(delta*b7)
	cov[3,3] = 24*b1/(I*b7)
	cov[3,4] = cov[4,3] = 24*b2/(I*b7)
	cov[4,4] = b10/(I*b7)
	end
	cov .*= sigma*sigma/Lambda_eff
	obs_dist = MvNormal(zeros(4),cov)

	local obs_duration, obs_depth, sigma_duration, sigma_depth
        isvalid = false
	while !isvalid
	  obs_vec = rand(obs_dist)
          obs_duration = duration + obs_vec[2]
          obs_depth = depth + obs_vec[3]
          if (obs_duration>0.0) && (obs_depth>0.0)
             isvalid = true
	  end
	end
     	    obs = TransitPlanetObs( period*(1.0+sigma_obs.period*randtn()), t0*(1.0+sigma_obs.period*randtn()), obs_depth,obs_duration)
        end
  	return obs, sigma_obs
end
