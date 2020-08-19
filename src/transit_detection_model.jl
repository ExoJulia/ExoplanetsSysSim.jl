## ExoplanetsSysSim/src/transit_detection_model.jl
## (c) 2015 Eric B. Ford

# Several functions below based on https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
# That follows the procedure outlined in Burke et al.(2015).
# However we don't currently interpolate the mesthreshold to the relevant duration
import SpecialFunctions.lgamma

function real_log_choose(m::Float64, n::Float64)::Float64
  lgamma(m+1)-lgamma(n+1)-lgamma(m-n+1.0)
end

function real_binom(k::Float64, BigM::Float64, f::Float64)::Float64
    F1 = real_log_choose(BigM,k)
    F2 = k*log(f)
    F3 = (BigM-k)*log(1.0-f)
    x = exp(F1+F2+F3)
    return x
end

"""
    kepler_window_function_binomial_model(exp_num_transits_no_gaps, duty_cycle; min_transits = 3.0)

Binomial window function model for use with Kepler taken from Burke et al. (2015).

# Arguments:
- `exp_num_transits_no_gaps::Float64`: Expected number of transits assuming no gaps in the observation
- `duty_cycle::Float64`: Kepler duty cycle (i.e. fraction of observation time with usable data)
- `min_transits::Float64 = 3.0`: Minimum number of observed transits required for consideration as planet candidate

# Returns:
Probability of detecting at least min_transits given the provided properties of the Kepler target star and planet.
"""
function kepler_window_function_binomial_model(exp_num_transits_no_gaps::Float64, duty_cycle::Float64; min_transits::Float64 = 3.0)::Float64
  if exp_num_transits_no_gaps < min_transits
     return 0.0
  else
     return max(1.0 - real_binom(min_transits,exp_num_transits_no_gaps,duty_cycle), 0.0)
  end
end

"""
    kepler_window_function_binomial_model(t, exp_num_transits_no_gaps, period, duration; min_transits = 3.0)

Binomial window function model for use with Kepler taken from Burke et al. (2015).

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `exp_num_transits_no_gaps::Float64`: Expected number of transits assuming no gaps in the observation
- `period::Float64`: Orbital period of planet (in days)
- `duration::Float64`: Transit duration of planet (in days)
- `min_transits::Float64 = 3.0`: Minimum number of observed transits required for consideration as planet candidate

# Returns:
Probability of detecting at least min_transits given the provided properties of the Kepler target star and planet.
"""
function kepler_window_function_binomial_model(t::KeplerTarget, exp_num_transits_no_gaps::Float64, period::Float64, duration::Float64; min_transits::Float64 = 3.0)::Float64
   kepler_window_function_binomial_model(exp_num_transits_no_gaps, t.duty_cycle, min_transits=min_transits)
end

"""
    kepler_window_function_dr25_model(t, exp_num_transits_no_gaps, period, duration)

Window function model for use with Kepler that were fit per target for DR25 (see KSCI-19101-002).

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `exp_num_transits_no_gaps::Float64`: Expected number of transits assuming no gaps in the observation
- `period::Float64`: Orbital period of planet (in days)
- `duration::Float64`: Transit duration of planet (in days)

# Returns:
Probability of detecting at least 3 transits given the provided properties of the Kepler target star and planet.
"""
function kepler_window_function_dr25_model(t::KeplerTarget, exp_num_transits_no_gaps::Float64, period::Float64, duration::Float64)::Float64
   ExoplanetsSysSim.WindowFunction.eval_window_function(t.window_function_id, Duration=duration, Period=period)
end

# WARNING: Hard-coded function variable for computational efficiency,
#          replace with appropriate window function model.
#kepler_window_function = kepler_window_function_binomial_model
kepler_window_function = kepler_window_function_dr25_model

function tess_window_function(t::TESSTarget, period::Float64, duration::Float64, num_transits_to_confirm::Int64 = 2)::Float64
    # if the window's too short to see one transit
    if period + duration > tess_sector_duration
        return num_transits_to_confirm == 1 ? tess_sector_duration / period : 0
    # if the window's guaranteed to have `num_transits_to_confirm` transits
    elseif floor((tess_sector_duration - period + duration) / (period + duration)) >= num_transits_to_confirm - 1
        return 1
    else
        # find P of t0 being small enough
        return (tess_sector_duration - (num_transits_to_confirm - 1) * period) / (period - duration)
    end
end

"""
    frac_depth_to_tps_depth(frac_depth)

Converts fractional transit depth to depth used by Kepler Transiting Planet Search (TPS) module. (See Christiansen et al. (2015))  
NOTE: This assumes a linear limbdarkening coefficient of 0.6

# Arguments:
- `frac_depth::Float64`: Fractional transit depth of planet at transit center.

# Returns:
Kepler TPS depth.
"""
function frac_depth_to_tps_depth(frac_depth::Float64)
    alp = 1.0874
    bet = 1.0187
    REALDEPTH2TPSSQUARE = 1.0  # WARNING: Waiting for this to be confirmed
    k = sqrt(frac_depth)
    tps_depth = min( (alp-bet*k) * frac_depth* REALDEPTH2TPSSQUARE, 1.0)   # NOTE: I added the max based on common sense
    return tps_depth::Float64
end

"""
    detection_efficiency_theory(mes, expected_num_transits; min_pdet_nonzero = 0.0)

Detection efficiency model assuming perfect theoretical 7.1-sigma error function response.

# Arguments: 
- `mes::Float64`: Estimated multiple event statistic (or signal-to-noise ratio) for planet
- `expected_num_transits::Float64`: Expected number of planet transits (assuming no observation gaps)
- `min_pdet_nonzero::Float64 = 0.0`: Minimum probability of detection (if transiting) to be treated as identifiable

# Results:
Probability of detection (if transiting) for planet
"""
function detection_efficiency_theory(mes::Float64, expected_num_transits::Float64; min_pdet_nonzero::Float64 = 0.0)
   muoffset =  0.0
   sig =  1.0
   mesthresh = 7.1
   mes *= 1.003
   if mes > (9.0 - mesthresh - muoffset)
      return 1.0
   else
      pdet = 0.5 + 0.5*erf((mes - mesthresh - muoffset) / sqrt(2.0*sig*sig))
      pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
      return pdet
   end
end

"""
    detection_efficiency_fressin2013(mes, expected_num_transits)

Detection efficiency model using the linear ramp from Fressin (2013).

# Arguments: 
- `mes::Float64`: Estimated multiple event statistic (or signal-to-noise ratio) for planet
- `expected_num_transits::Float64`: Expected number of planet transits (assuming no observation gaps)

# Results:
Probability of detection (if transiting) for planet
"""
function detection_efficiency_fressin2013(mes::Float64, expected_num_transits::Float64)
    mesmin =  6.0
    mesmax =  16.0
   if mes <= mesmin
      return 0.0
   elseif mes >= mesmax
      return 1.0
   else
     return (mes - mesmin) / (mesmax - mesmin)
  end
end

"""
    detection_efficiency_christiansen2015(mes, expected_num_transits; mes_threshold = 7.1, min_pdet_nonzero = 0.0)

Detection efficiency model using a Gamma function taken from Christiansen (2015).

# Arguments: 
- `mes::Float64`: Estimated multiple event statistic (or signal-to-noise ratio) for planet
- `expected_num_transits::Float64`: Expected number of planet transits (assuming no observation gaps)
- `mes_threshold::Float64 = 7.1`: Minimum multiple event statistic value for signal to have been considered a potential planet candidate signal
- `min_pdet_nonzero::Float64 = 0.0`: Minimum probability of detection (if transiting) to be treated as identifiable

# Results:
Probability of detection (if transiting) for planet
"""
function detection_efficiency_christiansen2015(mes::Float64, expected_num_transits::Float64; mes_threshold::Float64 = 7.1, min_pdet_nonzero::Float64 = 0.0)
    a =  4.65  # from code for detection_efficiency(...) at https://github.com/christopherburke/KeplerPORTs/blob/master/KeplerPORTs_utils.py
   # b =  1.05
   # a =  4.35 # from arxiv abstract. Informal testing showed it didn't matter
    b =  0.98
   mes *= 1.003
   usemes = max(0.0,mes - 7.1 - (mes_threshold - 7.1))
   pdet = cdf(Gamma(a,b), usemes)
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

"""
    detection_efficiency_dr25_simple(mes, expected_num_transits; min_pdet_nonzero = 0.0)

Detection efficiency model using a Gamma function fit over a FGK sample of DR25 targets, taken from KSCI-19110-001.

# Arguments: 
- `mes::Float64`: Estimated multiple event statistic (or signal-to-noise ratio) for planet
- `expected_num_transits::Float64`: Expected number of planet transits (assuming no observation gaps)
- `min_pdet_nonzero::Float64 = 0.0`: Minimum probability of detection (if transiting) to be treated as identifiable

# Results:
Probability of detection (if transiting) for planet
"""
function detection_efficiency_dr25_simple(mes::Float64, expected_num_transits::Float64; min_pdet_nonzero::Float64 = 0.0)::Float64
   a = 30.87  # from pg 16 of https://exoplanetarchive.ipac.caltech.edu/docs/KSCI-19110-001.pdf
   b = 0.271
   c = 0.940
   mes *= 1.003
   dist = Gamma(a,b)
   pdet::Float64 = c*cdf(dist, mes)::Float64
   pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
   return pdet
end

"""
    get_param_for_detection_and_vetting_efficiency_depending_on_num_transits(num_tr)

Determine appropriate parameters to use in Gamma function model for detection efficiency and vetting of Kepler planet candidates.  Fit using the simulated transit injection tests on Kepler DR25 targets.

# Arguments: 
- `num_tr::Integer`: (Expected) number of observed transits for planet

# Returns:
alpha, beta, and C, the three parameters to be used for the Gamma detection efficiency function 
"""
function get_param_for_detection_and_vetting_efficiency_depending_on_num_transits(num_tr::Integer)
    if num_tr <= 3
        return (33.3884, 0.264472, 0.699093)
    elseif num_tr <= 4
        return  (32.886, 0.269577, 0.768366)
    elseif num_tr <= 5
        return (31.5196, 0.282741, 0.833673)
    elseif num_tr <= 6
        return (30.9919, 0.286979, 0.859865)
    elseif num_tr <= 9
        return (30.1906, 0.294688, 0.875042)
    elseif num_tr <= 18
        return (31.6342, 0.279425, 0.886144)
    elseif num_tr <= 36
        return (32.6448, 0.268898, 0.889724)
    else
        return (27.8185, 0.32432, 0.945075)
    end
end

"""
    detection_and_vetting_efficiency_model_v1(mes, expected_num_transits; min_pdet_nonzero = 0.0)

Detection efficiency model using a Gamma function whose parameters were determined using a fit over the simulated transit injection tests of a sample of FGK DR25 targets (see Hsu et al. (2019)).
# WARNING: Combined detection and vetting efficiency model - do NOT include additional vetting efficiency

# Arguments:
- `mes::Float64`: Estimated multiple event statistic (or signal-to-noise ratio) for planet
- `expected_num_transits::Float64`: Expected number of planet transits (assuming no observation gaps)
- `min_pdet_nonzero::Float64 = 0.0`: Minimum probability of detection (if transiting) to be treated as identifiable

# Returns:
Probability of detection (if transiting) for planet
"""
function detection_and_vetting_efficiency_model_v1(mes::Float64, expected_num_transits::Float64; min_pdet_nonzero::Float64 = 0.0)::Float64
    mes *= 1.003
    num_transit_int = convert(Int64,floor(expected_num_transits))
    num_transit_int += rand() < expected_num_transits-num_transit_int ? 1 : 0
    a, b, c = get_param_for_detection_and_vetting_efficiency_depending_on_num_transits(num_transit_int)
    dist = Gamma(a,b)
    pdet::Float64 = c*cdf(dist, mes)::Float64
    pdet = pdet >= min_pdet_nonzero ? pdet : 0.0
    return pdet
end

# WARNING: Hardcoded choice of transit detection efficiency here for speed and so as to not have it hardcoded in multiple places
#detection_efficiency_model = detection_efficiency_christiansen2015
# detection_efficiency_model = detection_efficiency_dr25_simple
detection_efficiency_model = detection_and_vetting_efficiency_model_v1

"""
    vetting_efficiency_none(R_p, P)

Vetting efficiency model assuming all identified planet candidate signals are true positives

# Arguments:
- `R_p::Real`: Planet radius
- `P::Real`: Planet orbital period

# Returns:
1
"""
function vetting_efficiency_none(R_p::Real, P::Real)
    return 1.0
end

"""
    vetting_efficiency_dr25_mulders(R_p, P)

Vetting efficiency model from Mulders et al. (2018) using parameters fit for Gaia DR2 stellar properties applied to Kepler DR25 targets independent of reliability score (through e-mail communication).

# Arguments:
- `R_p::Real`: Planet radius (in R_sol)
- `P::Real`: Planet orbital period (in days)

# Returns:
Probability of planet with given properties to be vetted as a planet candidate.
"""
function vetting_efficiency_dr25_mulders(R_p::Real, P::Real)
    c = 0.93
    a_R = -0.03
    P_break = 205.
    a_P = 0.00
    b_P = -0.24
    pvet = c*(R_p/earth_radius)^a_R
    if P < P_break
        pvet *= (P/P_break)^a_P
    else
        pvet *= (P/P_break)^b_P
    end
    return pvet
end

"""
    vetting_efficiency_dr25_mulders_score_cut(R_p, P)

Vetting efficiency model from Mulders et al. (2018) using parameters fit for Gaia DR2 stellar properties applied to Kepler DR25 targets selected via reliability score (arXiv 1805.08211).

# Arguments:
- `R_p::Real`: Planet radius (in R_sol)
- `P::Real`: Planet orbital period (in days)

# Returns:
Probability of planet with given properties to be vetted as a planet candidate.
"""
function vetting_efficiency_dr25_mulders_score_cut(R_p::Real, P::Real)
    c = 0.63
    a_R = 0.19
    P_break = 53.
    a_P = -0.07
    b_P = -0.39
    pvet = c*(R_p/earth_radius)^a_R
    if P < P_break
        pvet *= (P/P_break)^a_P
    else
        pvet *= (P/P_break)^b_P
    end
    return pvet
end

# WARNING: Hardcoded choice of planet vetting efficiency here for speed and so as to not have it hardcoded in multiple places
vetting_efficiency = vetting_efficiency_none

"""
    interpolate_cdpp_to_duration_use_target_cdpp(t, duration)

Linearly interpolated CDPP using list of provided CDPP at different transit durations stored in Kepler target object (DEPRECATED).

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `duration::Float64`: Transit duration (in days)

# Results:
Estimated CDPP at given transit duration for Kepler target.
"""
function interpolate_cdpp_to_duration_use_target_cdpp(t::KeplerTarget, duration::Float64)::Float64
   duration_in_hours = duration *24.0
   dur_idx = searchsortedlast(cdpp_durations,duration_in_hours)   # cdpp_durations is defined in constants.jl
   if dur_idx <= 0
      cdpp = t.cdpp[1,1]
   elseif dur_idx==length(cdpp_durations) && (duration_in_hours >= cdpp_durations[end]) # Should be 15 cdpp_durations.
      cdpp = t.cdpp[length(cdpp_durations),1]
   else
      w = ((duration_in_hours)-cdpp_durations[dur_idx]) / (cdpp_durations[dur_idx+1]-cdpp_durations[dur_idx])
      cdpp = w*t.cdpp[dur_idx+1,1] + (1-w)*t.cdpp[dur_idx,1]
   end
   return cdpp
end

"""
    interpolate_cdpp_to_duration_lookup_cdpp(t, duration)

Linearly interpolated CDPP using list of provided CDPP at different transit durations stored in stellar catalog dataframe.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `duration::Float64`: Transit duration (in days)

# Results:
Estimated CDPP at given transit duration for Kepler target.
"""
function interpolate_cdpp_to_duration_lookup_cdpp(t::KeplerTarget, duration::Float64)::Float64
   duration_in_hours = duration *24.0
   dur_idx = searchsortedlast(cdpp_durations,duration_in_hours)   # cdpp_durations is defined in constants.jl
   get_cdpp(i::Integer) = 1.0e-6*sqrt(cdpp_durations[i]/24.0/kepler_LC_duration)*star_table(t,duration_symbols[i])

   if dur_idx <= 0
      cdpp = get_cdpp(1)
   elseif dur_idx==length(cdpp_durations) && (duration_in_hours >= cdpp_durations[end]) # Should be 15 cdpp_durations.
      cdpp = get_cdpp(length(cdpp_durations))
   else
      w = ((duration_in_hours)-cdpp_durations[dur_idx]) / (cdpp_durations[dur_idx+1]-cdpp_durations[dur_idx])
      cdpp = w*get_cdpp(dur_idx+1) + (1-w)*get_cdpp(dur_idx)
   end
   return cdpp
end

#interpolate_cdpp_to_duration = interpolate_cdpp_to_duration_use_target_cdpp
interpolate_cdpp_to_duration = interpolate_cdpp_to_duration_lookup_cdpp

"""
    calc_snr_if_transit_cdpp(t, depth, duration, cdpp, sim_param; num_transit = 1)

Calculate the expected multiple event statistic (signal-to-noise ratio) for planet around Kepler target star using CDPP.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `depth::Float64`: (Fractional) transit depth
- `duration::Float64`: Transit duration (in days)
- `cdpp::Float64`: CDPP for target star given transit duration
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns:
(Expected) multiple event statistic (signal-to-noise ratio)
"""
function calc_snr_if_transit_cdpp(t::KeplerTarget, depth::Float64, duration::Float64, cdpp::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  # depth_tps = frac_depth_to_tps_depth(depth)                  # TODO SCI:  WARNING: Hardcoded this conversion.  Remove once depth calculated using limb darkening model
  # snr = depth_tps*sqrt(num_transit*duration*kepler_LC_rate)/cdpp     # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
  snr = depth*sqrt(num_transit*duration*kepler_LC_rate)/cdpp     # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
end

function calc_snr_if_transit_cdpp(t::TESSTarget, depth::Float64, duration::Float64, noise::Float64, sim_param::SimParam; num_transit::Float64 = 1)
    # depth_tps = frac_depth_to_tps_depth(depth)                  # TODO SCI:  WARNING: Hardcoded this conversion.  Remove once depth calculated using limb darkening model
    # snr = depth_tps*sqrt(num_transit*duration*kepler_LC_rate)/cdpp     # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
    snr = depth*sqrt(num_transit*duration*tess_LC_rate)/noise     # WARNING: Assumes measurement uncertainties are uncorrelated & CDPP based on LC
  end

"""
    calc_snr_if_transit(t, depth, duration, osd, sim_param; num_transit = 1)

Calculate the expected multiple event statistic (signal-to-noise ratio) for planet around Kepler target star using 1-sigma depth function (OSD).
# NOTE: Assumes OSD functions have already been read in for all relevant Kepler targets.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `depth::Real`: (Fractional) transit depth
- `duration::Real`: Transit duration (in days)
- `osd::Real`: OSD for target star given transit duration and period
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Real = 1`: Expected number of transits observed

# Returns:
(Expected) multiple event statistic (signal-to-noise ratio)
"""
function calc_snr_if_transit(t::KeplerTarget, depth::Real, duration::Real, osd::Real, sim_param::SimParam; num_transit::Real = 1)
  # depth_tps = frac_depth_to_tps_depth(depth)                  # WARNING: Hardcoded this conversion
  # snr = depth_tps/osd*1.0e6 # osd is in ppm
  snr = depth/osd*1.0e6 # osd is in ppm
end

"""
    calc_snr_if_transit_central(t, s, p, sim_param)

Calculate the expected multiple event statistic (signal-to-noise ratio) for planet around Kepler target star using 1-sigma depth function (OSD) at transit center.
# NOTE: Assumes OSD functions have already been read in for all relevant Kepler targets.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns:
(Expected) multiple event statistic (signal-to-noise ratio) at transit center
"""
function calc_snr_if_transit_central(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  period = t.sys[s].orbit[p].P
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_eff_central(t,s,p)
  kepid = StellarTable.star_table(t.sys[s].star.id, :kepid)
  osd_duration_central = get_durations_searched_Kepler(period,duration_central)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
  osd_central = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration_central)
  if osd_duration_central > duration_central				#use a correcting factor if this duration is lower than the minimum searched for this period.
      osd_central = osd_central*osd_duration_central/duration_central
  end
  num_transit = calc_expected_num_transits(t,s,p,sim_param)
  calc_snr_if_transit(t,depth,duration_central,osd_central, sim_param,num_transit=num_transit)
end

"""
    calc_snr_if_transit_central_cdpp(t, s, p, sim_param)

Calculate the expected multiple event statistic (signal-to-noise ratio) for planet around Kepler target star using CDPP at transit center.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns:
(Expected) multiple event statistic (signal-to-noise ratio) at transit center
"""
function calc_snr_if_transit_central_cdpp(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_eff_central(t,s,p)
  cdpp = interpolate_cdpp_to_duration(t, duration_central)
  num_transit = calc_expected_num_transits(t,s,p,sim_param)
  calc_snr_if_transit_cdpp(t,depth,duration_central,cdpp, sim_param,num_transit=num_transit)
end

"""
    calc_prob_detect_if_transit(t, snr, period, duration, sim_param; num_transit = 1)
    calc_prob_detect_if_transit(t, depth, period, duration, osd, sim_param; num_transit = 1)

Calculate probability of detecting planet signal (if planet transits) using 1-sigma depth (OSD) function.
# NOTE: Assumes OSD functions already read in.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `snr::Float64`: Expected multiple event statistic (signal-to-noise ratio)
- `depth::Float64`: (Fractional) transit depth
- `period::Float64`: Orbital period (in days)
- `duration::Float64`: Transit duration (in days)
- `osd::Float64`: OSD for target star given transit duration and period
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns: 
Probability of detecting planet (if it transits)
"""
function calc_prob_detect_if_transit(t::KeplerTarget, snr::Float64, period::Float64, duration::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  min_pdet_nonzero = 1.0e-4                                                # TODO OPT: Consider raising threshold to prevent a plethora of planets that are very unlikely to be detected due to using 0.0 or other small value here
  wf = kepler_window_function(t, num_transit, period, duration)
  return wf*detection_efficiency_model(snr, num_transit, min_pdet_nonzero=min_pdet_nonzero)
end

function calc_prob_detect_if_transit(t::TESSTarget, snr::Array{Float64,1}, period::Float64, duration::Float64, sim_param::SimParam; num_transit::Float64 = 1)
    # WARNING: just returns 1, until prob_detect for TESS is properly defined
    return 1
end

function calc_prob_detect_if_transit(t::KeplerTarget, depth::Float64, period::Float64, duration::Float64, osd::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  snr = calc_snr_if_transit(t,depth,duration,osd, sim_param, num_transit=num_transit)
  return calc_prob_detect_if_transit(t, snr, period, duration, sim_param, num_transit=num_transit)
end

"""
    calc_prob_detect_if_transit_cdpp(t, depth, period, duration, cdpp, sim_param; num_transit = 1)

Calculate probability of detecting planet signal (if planet transits) using CDPP.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `depth::Float64`: (Fractional) transit depth
- `period::Float64`: Orbital period (in days)
- `duration::Float64`: Transit duration (in days)
- `cdpp::Float64`: CDPP for target star given transit duration
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns: 
Probability of detecting planet (if it transits)
"""
function calc_prob_detect_if_transit_cdpp(t::KeplerTarget, depth::Float64, period::Float64, duration::Float64, cdpp::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  snr = calc_snr_if_transit_cdpp(t,depth,duration,cdpp, sim_param, num_transit=num_transit)
  return calc_prob_detect_if_transit(t, snr, period, duration, sim_param, num_transit=num_transit)
end

"""
    calc_prob_detect_if_transit_central(t, s, p, sim_param)

Calculate probability of detecting planet signal (if planet transits) at transit center using 1-sigma depth (OSD) functions.
# NOTE: Assumes OSD functions already read in.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns: 
Probability of detecting planet (if it transits) at transit center
"""
function calc_prob_detect_if_transit_central(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  period = t.sys[s].orbit[p].P
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_eff_central(t,s,p)
  kepid = StellarTable.star_table(t.sys[s].star.id, :kepid)
  osd_duration_central = get_durations_searched_Kepler(period,duration_central)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
  osd_central = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration_central)
  if osd_duration_central > duration_central				#use a correcting factor if this duration is lower than the minimum searched for this period.
      osd_central = osd_central*osd_duration_central/duration_central
  end
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit(t,depth,period,duration_central,osd_central, sim_param, num_transit=ntr)
end

"""
    calc_prob_detect_if_transit_central_cdpp(t, s, p, sim_param)

Calculate probability of detecting planet signal (if planet transits) at transit center using CDPP

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns: 
Probability of detecting planet (if it transits) at transit center
"""
function calc_prob_detect_if_transit_central_cdpp(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  period = t.sys[s].orbit[p].P
  depth = calc_transit_depth(t,s,p)
  duration_central = calc_transit_duration_eff_central(t,s,p)
  cdpp = interpolate_cdpp_to_duration(t, duration_central)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit_cdpp(t,depth,period,duration_central,cdpp, sim_param, num_transit=ntr)
end

"""
    calc_prob_detect_if_transit_with_actual_b(t, s, p, sim_param)

Calculate probability of detecting planet signal (if planet transits) at transit center using 1-sigma depth (OSD) functions and the impact parameter of the provided orbit.
# NOTE: Assumes OSD functions already read in.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns: 
Probability of detecting planet (if it transits) given the impact parameter of the planet's orbit
"""
function calc_prob_detect_if_transit_with_actual_b(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  period = t.sys[s].orbit[p].P
  size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  depth = calc_transit_depth(t,s,p)
  duration = calc_transit_duration_eff(t,s,p)
  b = calc_impact_parameter(t.sys[s],p)
  snr_correction = calc_depth_correction_for_grazing_transit(b,size_ratio)
  depth *= snr_correction
  kepid = StellarTable.star_table(t.sys[s].star.id, :kepid)
  osd_duration = get_durations_searched_Kepler(period,duration)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
  osd = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration)
  if osd_duration > duration				#use a correcting factor if this duration is lower than the minimum searched for this period.
      osd = osd*osd_duration/duration
  end
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit(t,depth,period,duration,osd, sim_param, num_transit=ntr)
end

"""
    calc_prob_detect_if_transit_with_actual_b_cdpp(t, s, p, sim_param)

Calculate probability of detecting planet signal (if planet transits) at transit center using CDPP and the impact parameter of the provided orbit.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters

# Returns: 
Probability of detecting planet (if it transits) given the impact parameter of the planet's orbit
"""
function calc_prob_detect_if_transit_with_actual_b_cdpp(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  period = t.sys[s].orbit[p].P
  size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  depth = calc_transit_depth(t,s,p)
  duration = calc_transit_duration_eff(t,s,p)
  b = calc_impact_parameter(t.sys[s],p)
  snr_correction = calc_depth_correction_for_grazing_transit(b,size_ratio)
  depth *= snr_correction
  cdpp = interpolate_cdpp_to_duration(t, duration)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_prob_detect_if_transit_cdpp(t,depth,period,duration,cdpp, sim_param, num_transit=ntr)
end

"""
    calc_ave_prob_detect_if_transit_from_snr(t, snr_central, period, duration_central, size_ratio, osd_central, sim_param; num_transit = 1)

Calculate probability of detecting planet signal (if planet transits) using 1-sigma depth (OSD) function and averaged over impact parameters b~U[0,1).
# NOTE: Assumes OSD functions already read in.

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `snr_central::Float64`: Expected multiple event statistic (signal-to-noise ratio) at transit center
- `period::Float64`: Orbital period (in days)
- `duration_central::Float64`: Transit duration (in days) at transit center
- `size_ratio::Float64`: Ratio of planet-to-star radii
- `osd_central::Float64`: OSD for target star given transit duration at transit center and period
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns: 
Probability of detecting planet (if it transits) averaged over impact parameter
"""
function calc_ave_prob_detect_if_transit_from_snr(t::KeplerTarget, snr_central::Float64, period::Float64, duration_central::Float64, size_ratio::Float64, osd_central::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  min_pdet_nonzero = 1.0e-4
  wf = kepler_window_function(t, num_transit, period, duration_central)
  detection_efficiency_central = detection_efficiency_model(snr_central, num_transit, min_pdet_nonzero=min_pdet_nonzero)
  if wf*detection_efficiency_central <= min_pdet_nonzero
     return 0.
  end

  function calc_ave_prob_detect_if_transit_from_snr(t::TESSTarget, snr_central::Float64, period::Float64, duration_central::Float64, sim_param::SimParam, num_transit::Float64 = 1)
    # currently only uses info from one TESS sector at a time
    wf = tess_window_function(t, num_transit, period, duration_central)
    detection_efficiency_central = detection_efficiency_model(snr_central, num_transit, min_pdet_nonzero=min_pdet_nonzero)
  end
  # Breaking integral into two sections [0,1-b_boundary) and [1-b_boundary,1], so need at least 5 points to evaluate integral via trapezoid rule
  num_impact_param_low_b =  7                            # Number of points to evaluate integral over [0,1-b_boundary) via trapezoid rule
  num_impact_param_high_b = 5 # (size_ratio<=0.05) ? 5 : 11  # Number of points to evaluate integral over [1-b_boudnary,1) via trapezoid rule.  If using 2*size_ratio for bondary for small planets, then keep this odd, so one point lands on 1-size_ratio.
  @assert(num_impact_param_low_b >= 5)
  @assert(num_impact_param_high_b >= 3)
  num_impact_param = num_impact_param_low_b+num_impact_param_high_b-1 # One point is shared
  b_boundary = (size_ratio <= 0.15) ? 2*size_ratio : min(max(0.3,size_ratio),0.5)
  b = Array{Float64}(undef,num_impact_param)
  weight = Array{Float64}(undef,num_impact_param)
  b[1:num_impact_param_low_b] = range(0.0,stop=1-b_boundary,length=num_impact_param_low_b)
  b[num_impact_param_low_b:num_impact_param] .= range(1-b_boundary,stop=1.0,length=num_impact_param_high_b)
  weight[1:num_impact_param_low_b] .= (1-b_boundary)/(num_impact_param_low_b-1)  # Points for first integral
  weight[1] *= 0.5                        # Lower endpoint of first integral
  weight[num_impact_param_low_b] *= 0.5   # Upper endpoint of first integral
  weight[num_impact_param_low_b] += 0.5*(b_boundary)/(num_impact_param_high_b-1) # Also lower endpoint of second integral
  weight[(num_impact_param_low_b+1):num_impact_param] .= b_boundary/(num_impact_param_high_b-1)
  weight[num_impact_param] *= 0.5         # Upper endpoint of second integral
  @assert isapprox(sum(weight),1.0)

  function integrand(b::Float64)::Float64
     depth_factor = calc_depth_correction_for_grazing_transit(b,size_ratio)
     duration_factor = calc_transit_duration_eff_factor_for_impact_parameter_b(b,size_ratio)
     kepid = StellarTable.star_table(t.sys[1].star.id, :kepid)
     osd_duration = get_durations_searched_Kepler(period,duration_central*duration_factor)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
     osd = WindowFunction.interp_OSD_from_table(kepid, period, osd_duration)
     if osd_duration > duration_central*duration_factor				#use a correcting factor if this duration is lower than the minimum searched for this period.
	osd = osd*osd_duration/(duration_central*duration_factor)
     end
     snr_factor = depth_factor*(osd_central/osd)
     detection_efficiency_model(snr_central*snr_factor, num_transit, min_pdet_nonzero=min_pdet_nonzero)
  end

  ave_detection_efficiency = sum(weight .* map(integrand,b)::Vector{Float64} )

  return wf*ave_detection_efficiency
end

"""
    calc_ave_prob_detect_if_transit_from_snr_cdpp(t, snr_central, period, duration_central, size_ratio, osd_central, sim_param; num_transit = 1)

Calculate probability of detecting planet signal (if planet transits) using CDPP and averaged over impact parameters b~U[0,1).

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `snr_central::Float64`: Expected multiple event statistic (signal-to-noise ratio) at transit center
- `period::Float64`: Orbital period (in days)
- `duration_central::Float64`: Transit duration (in days) at transit center
- `size_ratio::Float64`: Ratio of planet-to-star radii
- `cdpp_central::Float64`: CDPP for target star given transit duration at transit center
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns: 
Probability of detecting planet (if it transits) averaged over impact parameter
"""
function calc_ave_prob_detect_if_transit_from_snr_cdpp(t::KeplerTarget, snr_central::Float64, period::Float64, duration_central::Float64, size_ratio::Float64, cdpp_central::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  min_pdet_nonzero = 1.0e-4
  wf = kepler_window_function(t, num_transit, period, duration_central)

  detection_efficiency_central = detection_efficiency_model(snr_central, num_transit, min_pdet_nonzero=min_pdet_nonzero)
  if wf*detection_efficiency_central <= min_pdet_nonzero
     return 0.
  end

  # Breaking integral into two sections [0,1-b_boundary) and [1-b_boundary,1], so need at least 5 points to evaluate integral via trapezoid rule
  num_impact_param_low_b =  7                            # Number of points to evaluate integral over [0,1-b_boundary) via trapezoid rule
  num_impact_param_high_b = 5 # (size_ratio<=0.05) ? 5 : 11  # Number of points to evaluate integral over [1-b_boudnary,1) via trapezoid rule.  If using 2*size_ratio for bondary for small planets, then keep this odd, so one point lands on 1-size_ratio.
  @assert(num_impact_param_low_b >= 5)
  @assert(num_impact_param_high_b >= 3)
  num_impact_param = num_impact_param_low_b+num_impact_param_high_b-1 # One point is shared
  b_boundary = (size_ratio <= 0.15) ? 2*size_ratio : min(max(0.3,size_ratio),0.5)
  b = Array{Float64}(undef,num_impact_param)
  weight = Array{Float64}(undef,num_impact_param)
  b[1:num_impact_param_low_b] = range(0.0,stop=1-b_boundary,length=num_impact_param_low_b)
  b[num_impact_param_low_b:num_impact_param] .= range(1-b_boundary,stop=1.0,length=num_impact_param_high_b)
  weight[1:num_impact_param_low_b] .= (1-b_boundary)/(num_impact_param_low_b-1)  # Points for first integral
  weight[1] *= 0.5                        # Lower endpoint of first integral
  weight[num_impact_param_low_b] *= 0.5   # Upper endpoint of first integral
  weight[num_impact_param_low_b] += 0.5*(b_boundary)/(num_impact_param_high_b-1) # Also lower endpoint of second integral
  weight[(num_impact_param_low_b+1):num_impact_param] .= b_boundary/(num_impact_param_high_b-1)
  weight[num_impact_param] *= 0.5         # Upper endpoint of second integral
  @assert isapprox(sum(weight),1.0)

  function integrand(b::Float64)::Float64
     depth_factor = calc_depth_correction_for_grazing_transit(b,size_ratio)
     duration_factor = calc_transit_duration_eff_factor_for_impact_parameter_b(b,size_ratio)
     cdpp = interpolate_cdpp_to_duration(t,duration_central*duration_factor)
     snr_factor = depth_factor*sqrt(duration_factor)*(cdpp_central/cdpp)
     detection_efficiency_model(snr_central*snr_factor, num_transit, min_pdet_nonzero=min_pdet_nonzero)
  end

  ave_detection_efficiency = sum(weight .* map(integrand,b)::Vector{Float64} )

  return wf*ave_detection_efficiency
end

"""
    calc_ave_prob_detect_if_transit_cdpp(t, depth, period, duration_central, size_ratio, sim_param; num_transit = 1)
    calc_ave_prob_detect_if_transit_cdpp(t, s, p, sim_param)

Calculate probability of detecting planet signal (if planet transits) using CDPP and averaged over impact parameters b~U[0,1).

# Arguments:
- `t::KeplerTarget`: Kepler target object
- `depth::Float64`: (Fractional) transit depth
- `period::Float64`: Orbital period (in days)
- `duration_central::Float64`: Transit duration (in days) at transit center
- `size_ratio::Float64`: Ratio of planet-to-star radii
- `s::Integer`: Star index within Kepler target
- `p::Integer`: Planet index within Kepler target
- `sim_param::SimParam`: Simulation parameters
- `num_transit::Float64 = 1`: Expected number of transits observed

# Returns: 
Probability of detecting planet (if it transits) averaged over impact parameter
"""
function calc_ave_prob_detect_if_transit_cdpp(t::KeplerTarget, depth::Float64, period::Float64, duration_central::Float64, size_ratio::Float64, sim_param::SimParam; num_transit::Float64 = 1)
  cdpp_central = interpolate_cdpp_to_duration(t, duration_central)
  snr_central = calc_snr_if_transit(t,depth,duration_central,cdpp_central, sim_param, num_transit=num_transit)
  return calc_ave_prob_detect_if_transit_from_snr_cdpp(t, snr_central, period, duration_central, size_ratio, cdpp_central, sim_param, num_transit=num_transit)
end

function calc_ave_prob_detect_if_transit_cdpp(t::KeplerTarget, s::Integer, p::Integer, sim_param::SimParam)
  size_ratio = t.sys[s].planet[p].radius/t.sys[s].star.radius
  depth = calc_transit_depth(t,s,p)
  period = t.sys[s].orbit[p].P
  duration_central = calc_transit_duration_eff_central(t,s,p)
  ntr = calc_expected_num_transits(t,s,p,sim_param)
  calc_ave_prob_detect_if_transit_cdpp(t,depth,period,duration_central, size_ratio, sim_param, num_transit=ntr)
end
