## ExoplanetsSysSim/src/target.jl
## (c) 2015 Eric B. Ford

#using Distributions

struct KeplerTarget
  #sys::PlanetarySystem   # Make array for planetary systems aroud multiple stars in one target?
  sys::Vector{PlanetarySystemAbstract}
  cdpp::Array{Float64,2} # fractional, not ppm; 2D to allow for multiple time scales, months, quarters or seasons/spacecraft rotation angles
                           # QUERY: Do we want a separate CDPP for SC?  Or will CDPP's in different months be for LC/SC depending on this variable?
                            # QUERY: Should this be moved to KeplerTargetObs?
                             # QUERY: Should we not add this to target and use the star id to lookup CDPP from the stellar table?
  contam::Float64     # QUERY: Do we want/need this, since we're able to generate multiple stars in a single target?
  data_span::Float64
  duty_cycle::Float64
  window_function_id::Int64 # Points to the id of the window function for this target
  #channel::Int64         # E.g., if we cared which Kepler channel the target fell on
  #has_sc::Vector{Bool}   # TODO OPT: Make Immutable Vector or BitArray for speed?  QUERY: Should this go in KeplerTargetObs?
  #                       # QUERY: Do we want a separate CDPP for SC?  Or will CDPP's in different months be for LC/SC depending on this variable?
  #ra::Float64           # E.g., if we cared about position on sky     QUERY:  Should we replace with galactic longitude and latitute?
  #dec::Floa64           #
end
num_planets(t::KeplerTarget) = sum( num_planets, t.sys)
flux(t::KeplerTarget) = sum(flux,t.sys)+t.contam

star_table(t::KeplerTarget, sym::Symbol) = StellarTable.star_table(t.sys[1].star.id,sym)


function draw_asymmetric_normal(mu::Real, sig_plus::Real, sig_minus::Real; rn = randn() )
  @assert sig_minus >= zero(sig_minus)
  mu + ( (stdn>=zero(stdn)) ? sig_plus*rn : sig_minus*rn )
end

function make_cdpp_array_empty(star_id::Integer)
  cdpp_arr = Array{Float64,2}(undef,0,0)
end

function make_cdpp_array(star_id::Integer)
  star_table(id::Integer,sym::Symbol) = StellarTable.star_table(id,sym)::Float64
  cdpp_arr = (1.0e-6*sqrt(1.0/24.0/LC_duration)) .* Float64[star_table(star_id, :rrmscdpp01p5)*sqrt(1.5), star_table(star_id, :rrmscdpp02p0)*sqrt(2.), star_table(star_id,:rrmscdpp02p5)*sqrt(2.5), star_table(star_id,:rrmscdpp03p0)*sqrt(3.), star_table(star_id,:rrmscdpp03p5)*sqrt(3.5), star_table(star_id,:rrmscdpp04p5)*sqrt(4.5), star_table(star_id,:rrmscdpp05p0)*sqrt(5.), star_table(star_id,:rrmscdpp06p0)*sqrt(6.), star_table(star_id,:rrmscdpp07p5)*sqrt(7.5), star_table(star_id,:rrmscdpp09p0)*sqrt(9.), star_table(star_id,:rrmscdpp10p5)*sqrt(10.5), star_table(star_id,:rrmscdpp12p0)*sqrt(12.), star_table(star_id,:rrmscdpp12p5)*sqrt(12.5), star_table(star_id,:rrmscdpp15p0)*sqrt(15.)]
end

function generate_kepler_target_from_table(sim_param::SimParam)
  # generate_star = get_function(sim_param,"generate_star")
  generate_planetary_system = get_function(sim_param,"generate_planetary_system")
  max_draws_star_properties = 20
  min_star_radius = 0.5
  min_star_mass = 0.5
  max_star_radius = 2.0
  max_star_mass = 2.0
  max_star_density = 1000.0
  use_star_table_sigmas = true
  min_frac_rad_sigma = 0.06
  max_star_id = StellarTable.num_usable_in_star_table()

  star_table(id::Integer,sym::Symbol) = StellarTable.star_table(id,sym)

  @assert(1<=max_star_id)
  star_id = rand(1:max_star_id)
  mass = 0.0
  dens = 0.0
  radius = 0.0
    #if use_star_table_sigmas
    if get(sim_param,"use_star_table_sigmas",false)
        attmpt_num = 0
        while (!(min_star_radius<radius<max_star_radius)) || (!(min_star_mass<mass<max_star_mass))# || (!(0.0<dens<max_star_density))
            if attmpt_num >= max_draws_star_properties
                star_id = rand(1:max_star_id)
                attmpt_num = 0
            end
            rad_errp = max(star_table(star_id,:radius_err1), min_frac_rad_sigma*star_table(star_id,:radius))
            rad_errn = max(abs(star_table(star_id,:radius_err2)), min_frac_rad_sigma*star_table(star_id,:radius))
            rn = randn()
            radius = draw_asymmetric_normal( star_table(star_id,:radius), rad_errp,  rad_errn, rn=rn)
            mass = draw_asymmetric_normal( star_table(star_id,:mass), star_table(star_id,:mass_err1), abs(star_table(star_id,:mass_err2)), rn=rn )
            #dens = draw_asymmetric_normal( star_table(star_id,:dens), star_table(star_id,:dens_err1), abs(star_table(star_id,:dens_err2)) )
            attmpt_num += 1
        end

        # # ZAMS mass-radius relation taken from 15.1.1 of Allen's Astrophysical Quantities (2002)
        # if radius > 1.227
        #     mass = 10^((log10(radius)-0.011)/0.64)
        # else
        #     mass = 10^((log10(radius)+0.02)/0.917)
        # end

        dens   = (mass*sun_mass_in_kg_IAU2010*1000.)/(4//3*pi*(radius*sun_radius_in_m_IAU2015*100.)^3)  # Self-consistent density (gm/cm^3)
  else
    radius = star_table(star_id,:radius)
    mass   = star_table(star_id,:mass)
    #dens   = star_table(star_id,:dens)
    dens   = (mass*sun_mass_in_kg_IAU2010*1000.)/(4//3*pi*(radius*sun_radius_in_m_IAU2015*100.)^3)  # Self-consistent density (gm/cm^3)
  end
  ld = LimbDarkeningParam4thOrder(star_table(star_id,:limbdark_coeff1), star_table(star_id,:limbdark_coeff2), star_table(star_id,:limbdark_coeff3), star_table(star_id,:limbdark_coeff4) )
  star = SingleStar(radius,mass,1.0,ld,star_id)     # TODO SCI: Allow for blends, binaries, etc.
  #cdpp_arr = make_cdpp_array(star_id)
  cdpp_arr = make_cdpp_array_empty(star_id) # TODO: WARNING: Now leaving this field empty out and looking up each time
  contam = star_table(star_id, :contam)
  data_span = star_table(star_id, :dataspan)
  duty_cycle = star_table(star_id, :dutycycle)
  if StellarTable.star_table_has_key(:wf_id)
     wf_id = star_table(star_id,:wf_id)
  else
     wf_id = WindowFunction.get_window_function_id(star_table(star_id,:kepid))
  end
  # ch = rand(DiscreteUniform(1,84))              # Removed channel in favor of window function id
  ps = generate_planetary_system(star, sim_param)
  return KeplerTarget([ps],repeat(cdpp_arr, outer=[1,1]),contam,data_span,duty_cycle,wf_id)
end

function generate_kepler_target_simple(sim_param::SimParam)
   generate_star = get_function(sim_param,"generate_star")
   generate_planetary_system = get_function(sim_param,"generate_planetary_system")
   star::StarAbstract = generate_star(sim_param)

   mean_log_cdpp = 4.9759601617565465   # mean frmo star table
   stddev_log_cdpp = 0.6704860437536709 # std dev from star table
  rrmscdpp_5hr = exp(mean_log_cdpp+stddev_log_cdpp*randn())
  cdpp_5hr = 1.0e-6 * rrmscdpp_5hr * sqrt(5.0/24.0 / LC_duration )
  contam = 0.0 # rand(LogNormal(1.0e-3,1.0))   # TODO SCI: Come up with better description of Kepler targets, maybe draw from real contaminations
  wf_id = 0
  # ch = rand(DiscreteUniform(1,84))              # Removed channel in favor of window function id
  ps = generate_planetary_system(star, sim_param)
  return KeplerTarget(PlanetarySystemAbstract[ps],fill(cdpp_5hr,num_cdpp_timescales,num_quarters),contam,mission_data_span,mission_duty_cycle,0)
end

function test_target(sim_param::SimParam)
  generate_kepler_target_simple(sim_param)
  StellarTable.setup_star_table(sim_param)
  generate_kepler_target_from_table(sim_param)
end
