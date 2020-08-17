## ExoplanetsSysSim/src/mission_constants.jl
## (c) 2015 Eric B. Ford

#module KeplerConstants
 const num_channels = 84
 const num_modules = 42
 const num_quarters = 17              # QUERY:  I'm favoring leaving out quarter 0, since that was engineering data.  Agree?
 const cdpp_durations = [1.5,2.,2.5,3.,3.5,4.5,5.,6.,7.5,9.,10.5,12.,12.5,15.]
 const duration_symbols = [:rrmscdpp01p5, :rrmscdpp02p0,:rrmscdpp02p5,:rrmscdpp03p0,:rrmscdpp03p5,:rrmscdpp04p5,:rrmscdpp05p0,:rrmscdpp06p0,:rrmscdpp07p5,:rrmscdpp09p0,:rrmscdpp10p5,:rrmscdpp12p0,:rrmscdpp12p5,:rrmscdpp15p0 ]
 const num_cdpp_timescales = 14
 @assert num_cdpp_timescales == length(cdpp_durations) == length(duration_symbols)
 const mission_data_span = 1459.789   # maximum(ExoplanetsSysSim.StellarTable.df[:dataspan])
 const mission_duty_cycle = 0.8751    # median(ExoplanetsSysSim.StellarTable.df[:dutycycle])

 const kepler_exp_time_internal  =  6.019802903/(24*60*60)    # https://archive.stsci.edu/kepler/manuals/archive_manual.pdf
 const kepler_read_time_internal = 0.5189485261/(24*60*60)    # https://archive.stsci.edu/kepler/manuals/archive_manual.pdf
 const num_exposures_per_LC = 270
 const num_exposures_per_SC = 9
 const LC_integration_time = kepler_exp_time_internal*num_exposures_per_LC
 const SC_integration_time = kepler_exp_time_internal*num_exposures_per_SC
 const LC_read_time = kepler_read_time_internal*num_exposures_per_LC
 const SC_read_time = kepler_read_time_internal*num_exposures_per_SC
 const LC_duration = LC_integration_time +  LC_read_time
 const SC_duration = SC_integration_time +  SC_read_time
 const LC_rate = 1.0/LC_duration
 const SC_rate = 1.0/SC_duration

 const tess_sector_duration = 27.4 # days

# export num_channels,num_modules,num_quarters,num_cdpp_timescales,mission_data_span,mission_duty_cycle
# export kepler_exp_time_internal,kepler_read_time_internal,num_exposures_per_LC,num_exposures_per_SC
# export LC_integration_time,SC_integration_time,LC_read_time,SC_read_time,LC_duration,SC_duration,LC_rate,SC_rate
#end

#using KeplerConstants

# Standard conversion factors on which unit system is based
const global AU_in_m_IAU2012 = 149597870700.0
 const global G_in_mks_IAU2015 = 6.67384e-11
 const global G_mass_sun_in_mks = 1.3271244e20
 const global G_mass_earth_in_mks = 3.986004e14
 const global sun_radius_in_m_IAU2015 = 6.9566e8
 const global earth_radius_eq_in_m_IAU2015 = 6.3781e6
 const global sun_mass_in_kg_IAU2010 = 1.988547e30

# Constants used by this code
const global sun_mass = 1.0
 const global earth_mass = G_mass_earth_in_mks/G_mass_sun_in_mks  # about 3.0024584e-6
 const global earth_radius = earth_radius_eq_in_m_IAU2015 / sun_radius_in_m_IAU2015 # about 0.0091705248
 const global rsol_in_au = sun_radius_in_m_IAU2015 / AU_in_m_IAU2012  # about 0.00464913034
 const global sec_in_day = 24*60*60
 const global grav_const = G_in_mks_IAU2015 * sec_in_day^2 * sun_mass_in_kg_IAU2010 / AU_in_m_IAU2012^3 # about 2.9591220363e-4 in AU^3/(day^2 Msol)
 const global day_in_year = 365.25
