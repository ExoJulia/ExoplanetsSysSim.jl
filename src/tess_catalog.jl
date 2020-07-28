## ExoplanetsSysSim/src/tess_catalog.jl
## (c) 2020 Aditya Sengupta?

using DataFrames
using CSV
using FileIO

mutable struct TESSPhysicalCatalog
  target::Array{TessTarget,1}
end

mutable struct TESSObsCatalog
  target::Array{TESSTargetObs,1}
end
TESSObsCatalog() = TESSObsCatalog(TESSTargetObs[])

"""
    generate_tess_physical_catalog(sim_param)

Wrapper function to create catalog of simulated TESS targets.

# Arguments:
- `sim_param::SimParam`: Simulation parameter object; must have the following parameters set:
  - num_targets_sim_pass_one = Number of TESS targets in simulated catalog
  - generate_tess_target = Function which generates TESS targets
  - (Optional) stellar_catalog = Stellar catalog filename
  - (Optional) star_table_setup = Function that loads stellar catalog into DataFrame

# Returns:
TESS physical catalog object containing all simulated TESS targets and associated planetary systems.
"""
function generate_tess_physical_catalog(sim_param::SimParam)
   if haskey(sim_param,"stellar_catalog")
      star_tab_func = get_function(sim_param, "star_table_setup")
      star_tab_func(sim_param)
   end
   num_sys = get_int(sim_param,"num_targets_sim_pass_one")
   generate_tess_target = get_function(sim_param,"generate_tess_target")
   target_list = Array{TESSTarget}(undef,num_sys)
   map!(x->generate_tess_target(sim_param), target_list, 1:num_sys )
   return TESSPhysicalCatalog(target_list)
end

"""
    observe_tess_targets_sky_avg(input, sim_param)

Wrapper function to create catalog of simulated observations of TESS targets using the sky averaging observation scheme (i.e. each planet's detection probability is the average detection probability over all view-angles).

# Arguments:
- `input::TESSPhysicalCatalog`: Catalog object of simulated TESS targets and associated planetary systems to be observed
- `sim_param::SimParam`: Simulation parameter object; requires the following simulation parameters to be set:
  - calc_target_obs_sky_ave: Function name for sky averaging simulated observations

# Returns:
TESS observations catalog object containing all properties observed from the TESS targets and associated planetary systems that were detected during the simulation.
"""
function observe_tess_targets_sky_avg(input::TESSPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  return observe_tess_targets(calc_target_obs, input, sim_param)
end

"""
    observe_tess_targets_single_obs(input, sim_param)

Wrapper function to create catalog of simulated observations of TESS targets using the single observer observation scheme (i.e. each planet's detection probability is the detection probability from the Earth).

# Arguments:
- `input::TESSPhysicalCatalog`: Catalog object of simulated TESS targets and associated planetary systems to be observed
- `sim_param::SimParam`: Simulation parameter object; requires the following simulation parameters to be set:
  - calc_target_obs_single_obs: Function name for single observer simulated observations

# Returns:
TESS observations catalog object containing all properties observed from the TESS targets and associated planetary systems that were detected during the simulation.
"""
function observe_tess_targets_single_obs(input::TESSPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  return observe_tess_targets(calc_target_obs, input, sim_param)
end

"""
    observe_tess_targets(calc_target_obs, input, sim_param)

Wrapper function to create catalog of simulated observations of TESS targets.

# Arguments:
- `calc_target_obs::Function`: Function to use in simulating observations of TESS targets (sky averaging vs. single observer schemes)
- `input::TESSPhysicalCatalog`: Catalog object of simulated TESS targets and associated planetary systems to be observed
- `sim_param::SimParam`: Simulation parameter object

# Returns:
TESS observations catalog object containing all properties observed from the TESS targets and associated planetary systems that were detected during the simulation.
"""
function observe_tess_targets(calc_target_obs::Function, input::TESSPhysicalCatalog, sim_param::SimParam )
  #calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  #calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  output = TESSObsCatalog()
  if haskey(sim_param,"mem_tess_target_obs")
     output.target = get(sim_param,"mem_tess_target_obs",Array{TESSTargetObs}(0) )
  end
  num_targets_sim_pass_one = get_int(sim_param,"num_targets_sim_pass_one")
  if length(output.target) < num_targets_sim_pass_one
     output.target = Array{TESSTargetObs}(undef,num_targets_sim_pass_one)
  end
  map!(x::TESSTarget->calc_target_obs(x,sim_param)::TESSTargetObs, output.target, input.target)
  resize!(output.target,length(input.target))
  return output
end

# Test if this planetary system has at least one planet that transits (assuming a single observer)
function select_targets_one_obs(ps::PlanetarySystemAbstract)
 for pl in 1:length(ps.orbit)
   ecc::Float64 = ps.orbit[pl].ecc
   incl::Float64 = ps.orbit[pl].incl
   a::Float64 = semimajor_axis(ps,pl)
   Rstar::Float64 = rsol_in_au*ps.star.radius
   if (Rstar > (a*(1-ecc)*(1+ecc))/(1+ecc*sin(ps.orbit[pl].omega))*cos(incl))
     return true
   end
 end
 return false
end
#=
function select_targets_one_obs(ps::PlanetarySystemAbstract)
  for pl in 1:length(ps.orbit)
    if does_planet_transit(ps,pl)
       return true
    end
  end
  return false
end
=#

# Remove undetected planets from physical catalog
# TODO: OPT: Maybe create array of bools for which planets to keep, rather than splicing out non-detections?
function generate_obs_targets(cat_phys::TESSPhysicalCatalog, sim_param::SimParam )
  for t in 1:length(cat_phys.target)
    kep_targ = cat_phys.target[t]
    for ps in 1:length(cat_phys.target[t].sys)
      sys = kep_targ.sys[ps]
      for pl in length(sys.orbit):-1:1    # Going in reverse since removing planets from end of list first is cheaper than starting at beginning
        ecc::Float64 = sys.orbit[pl].ecc
	incl::Float64 = sys.orbit[pl].incl
   	a::Float64 = semimajor_axis(sys,pl)
   	Rstar::Float64 = rsol_in_au*sys.star.radius

        does_it_transit = does_planet_transit(sys, pl)
        pdet_if_tr = does_it_transit ? calc_prob_detect_if_transit_with_actual_b(kep_targ, ps, pl, sim_param) : 0.
        if !does_it_transit || (rand()>pdet_if_tr)
    	  splice!(cat_phys.target[t].sys[ps].orbit, pl)
	  splice!(cat_phys.target[t].sys[ps].planet, pl)
     	end
      end
    end
  end
  return cat_phys
end


# The following function is primarily left for debugging.
# Create a catalog of observations of simulated TESS targets.
function simulated_read_tess_observations(sim_param::SimParam )
   println("# WARNING: Using simulated_read_tess_observations.")
   # if haskey(sim_param,"stellar_catalog")
   #    star_tab_func = get_function(sim_param, "star_table_setup")
   #    star_tab_func(sim_param)
   # end
   num_sys = get_int(sim_param,"num_tess_targets")
   generate_tess_target = get_function(sim_param,"generate_tess_target")
   target_list = Array{TESSTarget}(undef,num_sys)
   map!(x->generate_tess_target(sim_param), target_list, 1:num_sys )

   cat_phys_cut = generate_obs_targets(TESSPhysicalCatalog(target_list), sim_param)
   calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
   output = TESSObsCatalog()
   output.target = map(x::TESSTarget->calc_target_obs(x,sim_param)::TESSTargetObs, cat_phys_cut.target)
   return output
end

"""
    read_toi_catalog(sim_param, force_reread=false)

Wrapper function to read TESS Object of Interest (TOI) catalog given SimParam

# Arguments:
- `sim_param::SimParam`: Simulation parameter object; this function uses the following parameters from the SimParam object:
  - toi_catalog: String filename of TESS Object of Interest catalog
- `force_reread::Bool`: Should the file be read in even if a DataFrame of the TOIs already exists in workspace?

# Returns:
- Dataframe of TOI objects (and their respective properties).
- Vector of booleans indicating which TOIs  were designated as planet candidates by the TESS pipeline and have a valid observed radius ratio and period (necessary for detection probability calculation).
"""
function read_toi_catalog(sim_param::SimParam, force_reread::Bool = false)
    filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)),"..", "data", convert(String,get(sim_param,"toi_catalog")) ) )
    return read_toi_catalog(filename, force_reread)
end

"""
    read_toi_catalog(filename, force_reread=false)

Function to read TESS Object of Interest (TOI) catalog given filename string.

# Arguments:
- `filename::String`: String filename of TESS Object of Interest catalog
- `force_reread::Bool`: Should the file be read in even if a DataFrame of the TOIs already exists in workspace?

# Returns:
- Dataframe of TOI objects (and their respective properties).
- Vector of booleans indicating which TOIs were designated as planet candidates by the TESS pipeline and have a valid observed radius ratio and period (necessary for detection probability calculation).
"""
function read_toi_catalog(filename::String, force_reread::Bool = false)
    local df, usable

    if occursin(r".jld2$",filename) && !force_reread
        try
            data = load(filename)
            df = data["toi_catalog"]
            usable = data["toi_catalog_usable"]
            Core.typeassert(df,DataFrame)
            Core.typeassert(usable,Array{Int64,1})
        catch
            error(string("# Failed to read toi catalog >",filename,"< in jld2 format."))
        end
    else
       try
            df = CSV.read(filename,comment="#")

            # Choose which TOIs to keep
            #is_cand = (csv_data[!,:,toi_disposition_idx] .== "CONFIRMED") | (csv_data[!,:,toi_disposition_idx] .== "CANDIDATE")
            is_cand = df[!,:toi_pdisposition] .== "CANDIDATE"
            has_radius = .!ismissing.(df[!,:toi_ror])
            has_period = .!(ismissing.(df[!,:toi_period]) .| ismissing.(df[!,:toi_period_err1]) .| ismissing.(df[!,:toi_period_err2]))

            is_usable = .&(is_cand, has_radius, has_period)
            usable = findall(is_usable)
           #  symbols_to_keep = [:ticid, :kepoi_name, :toi_pdisposition, :toi_score, :toi_ror, :toi_period, :toi_period_err1, :toi_period_err2, :toi_time0bk, :toi_time0bk_err1, :toi_time0bk_err2, :toi_depth, :toi_depth_err1, :toi_depth_err2, :toi_duration, :toi_duration_err1, :toi_duration_err2]
           # df = df[usable, symbols_to_keep]
           # tmp_df = DataFrame()
           # for col in names(df)
           #     tmp_df[col] = collect(skipmissing(df[col]))
           # end
           # df = tmp_df
           # usable = collect(1:length(df[!,:ticid]))
        catch
            error(string("# Failed to read toi catalog >",filename,"< in ascii format."))
        end
    end
    return df, usable
end

"""
    setup_actual_planet_candidate_catalog(df_star, df_toi, usable_toi, sim_param)

Create (true) catalog of TESS observations of TESS targets

# Arguments:
- `df_star::DataFrame`: DataFrame containing all TESS target stars in catalog
   NOTE: df_star is assumed to have fields ticid, mass and radius for all targets in the survey)
- `df_toi::DataFrame`: DataFrame containing all TESS Object of Interests (TOIs)
- `usable_toi::Array{Integer}`: Array of TOI dataframe row indices corresponding to TOIs to use
- `sim_param::SimParam`: Simulation parameter object

# Returns:
- TESS observations catalog containing TESS targets and associated TOIs (to be used as true catalog in comparison with simulated observations)
"""
function setup_actual_planet_candidate_catalog(df_star::DataFrame, df_toi::DataFrame, usable_toi::Array{Int64}, sim_param::SimParam)
    local target_obs, num_pl
    df_toi = df_toi[usable_toi,:]

    # Deprecated code to take a list of KepIDs and TOI names to pre-select a subset of TOIs
    # if haskey(sim_param, "toi_subset_csv")
    #     toi_subset = fill(false, length(df_toi[!,:ticid]))

    #     subset_df = readtable(convert(String,get(sim_param,"toi_subset_csv", "christiansen_kov.csv")), header=true, separator=' ')

    #     for n in 1:length(subset_df[!,1])
    #         subset_colnum = 1
    #         subset_entry = findall(x->x==subset_df[n,1], df_toi[names(subset_df)[1]])
    #         # println("Initial cut: ", subset_entry)
    #         while (length(subset_entry) > 1) & (subset_colnum < length(names(subset_df)))
    #             subset_colnum += 1

    #             subsubset = findall(x->round(x*10.)==round(subset_df[n,subset_colnum]*10.), df_toi[subset_entry,names(subset_df)[subset_colnum]])
    #             # println("Extra cut: ", subset_df[n,subset_colnum], " / ", df_toi[subset_entry,col_idx], " = ", subsubset)
    #             subset_entry = subset_entry[subsubset]
    #         end
    #         if length(subset_entry) > 1
    #             cand_sub = findall(x->x == "CANDIDATE",df_toi[subset_entry,:toi_pdisposition])
    #             subset_entry = subset_entry[cand_sub]
    #             if length(subset_entry) > 1
    #                 println("# Multiple planets found in final cut: ", subset_df[n,1])
    #             end
    #         end
    #         if length(subset_entry) < 1
    #             println("# No planets found in final cut: ", subset_df[n,:])
    #         end
    #         toi_subset[subset_entry] = true
    #     end
    #     df_toi = df_toi[findall(toi_subset),:]
    #     tot_plan = count(x->x, toi_subset)
    # end

    output = TESSObsCatalog()
    sort!(df_star, (:ticid))
    df_obs = join(df_star, df_toi, on = :ticid)
    #df_obs = sort!(df_obs, cols=(:ticid))
    df_obs = sort!(df_obs, (:ticid))

    # if haskey(sim_param, "toi_subset_csv")
    #     tot_plan -= length(df_obs[!,:kepoi_name])
    #     println("# Number of planet candidates in subset file with no matching star in table: ", tot_plan)
    # end

    # Add each TOI planet candidate to TESS target object
    plid = 0
    for i in 1:length(df_obs[!,:kepoi_name])
        if plid == 0
            plid = 1
            while i+plid < length(df_obs[!,:kepoi_name]) && df_obs[i+plid,:ticid] == df_obs[i,:ticid]
                plid += 1
            end
            num_pl = plid
            target_obs = TESSTargetObs(num_pl)
	        #target_obs.star = ExoplanetsSysSim.StarObs(df_obs[i,:radius],df_obs[i,:mass],findfirst(df_star[!,:ticid], df_obs[i,:ticid]))
            star_idx = searchsortedfirst(df_star[!,:ticid],df_obs[i,:ticid])
            if star_idx > length(df_star[!,:ticid])
                @warn "# Couldn't find ticid " * df_star[i,:ticid] * " in df_obs."
                star_idx = rand(1:length(df_star[!,:ticid]))
            end
            target_obs.star = ExoplanetsSysSim.StarObs(df_obs[i,:radius],df_obs[i,:mass],star_idx)

        end

        target_obs.obs[plid] = ExoplanetsSysSim.TransitPlanetObs(df_obs[i,:toi_period],df_obs[i,:toi_time0bk],df_obs[i,:toi_depth]/1.0e6,df_obs[i,:toi_duration])
        target_obs.sigma[plid] = ExoplanetsSysSim.TransitPlanetObs((abs(df_obs[i,:toi_period_err1])+abs(df_obs[i,:toi_period_err2]))/2,(abs(df_obs[i,:toi_time0bk_err1])+abs(df_obs[i,:toi_time0bk_err2]))/2,(abs(df_obs[i,:toi_depth_err1]/1.0e6)+abs(df_obs[i,:toi_depth_err2]/1.0e6))/2,(abs(df_obs[i,:toi_duration_err1])+abs(df_obs[i,:toi_duration_err2]))/2)
	#target_obs.prob_detect = ExoplanetsSysSim.SimulatedSystemDetectionProbs{OneObserver}( ones(num_pl), ones(num_pl,num_pl), ones(num_pl), fill(Array{Int64}(undef,0), 1) )  # Made line below to simplify calling
        target_obs.prob_detect = ExoplanetsSysSim.OneObserverSystemDetectionProbs(num_pl)
        plid -= 1
        if plid == 0
            push!(output.target,target_obs)
        end
    end
    return output
end

# Two functions below were just for debugging purposes
# Calculate SNR of every planet in simulated catalog
function calc_snr_list(cat::TESSPhysicalCatalog, sim_param::SimParam)
  snrlist = Array{Float64}(undef,0)
  for t in 1:length(cat.target)
    for p in 1:length(cat.target[t].sys[1].planet)
      snr = calc_snr_if_transit(cat.target[t],1,p,sim_param)
      if snr>0.0
        push!(snrlist,snr)
      end
    end
  end
  snrlist[findall(x->x>7.1,snrlist)]
end

# Calculate detection probability (assuming planet transits) for every planet in simulated catalog
function calc_prob_detect_list(cat::TESSPhysicalCatalog, sim_param::SimParam)
  pdetectlist = Array{Float64}(undef,0)
  for t in 1:length(cat.target)
    for p in 1:length(cat.target[t].sys[1].planet)
      #pdet = calc_prob_detect_if_transit(cat.target[t],1,p,sim_param)
      pdet = calc_prob_detect_if_transit_with_actual_b(cat.target[t],1,p,sim_param)

      if pdet>0.0
        push!(pdetectlist,pdet)
      end
    end
  end
  idx = findall(x->x>0.0,pdetectlist)
  pdetectlist[idx]
end

function test_catalog_constructors(sim_param::SimParam)
  cat_phys = generate_tess_physical_catalog(sim_param)::TESSPhysicalCatalog
  id = findfirst( x->num_planets(x)>=1 , cat_phys.target)   # fast forward to first target that has some planets
  @assert(length(id)>=1)
  semimajor_axis(cat_phys.target[id].sys[1],1)
  pdetlist = calc_prob_detect_list(cat_phys,sim_param)
  calc_target_obs_single_obs(cat_phys.target[id],sim_param)
  calc_target_obs_sky_ave(cat_phys.target[id],sim_param)
  @assert( length(cat_phys.target[id].sys[1].planet)  == num_planets(cat_phys.target[id]) )
  cat_obs = simulated_read_tess_observations(sim_param)
  return (cat_phys, cat_obs)
end
