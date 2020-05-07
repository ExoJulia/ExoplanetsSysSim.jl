## ExoplanetsSysSim/src/kepler_catalog.jl
## (c) 2015 Eric B. Ford

#using ExoplanetsSysSim
using DataFrames
#using DataArrays
using CSV
#using JLD
#using JLD2
using FileIO

#if VERSION >= v"0.5-"
#  import Compat: UTF8String, ASCIIString
#end

mutable struct KeplerPhysicalCatalog
  target::Array{KeplerTarget,1}
end
#KeplerPhysicalCatalog() = KeplerPhysicalCatalog([])

mutable struct KeplerObsCatalog
  target::Array{KeplerTargetObs,1}
end
KeplerObsCatalog() = KeplerObsCatalog(KeplerTargetObs[])

function generate_kepler_physical_catalog(sim_param::SimParam)
   if haskey(sim_param,"stellar_catalog")
      star_tab_func = get_function(sim_param, "star_table_setup")
      star_tab_func(sim_param)
   end
   num_sys = get_int(sim_param,"num_targets_sim_pass_one")
   generate_kepler_target = get_function(sim_param,"generate_kepler_target")
   target_list = Array{KeplerTarget}(undef,num_sys)
   map!(x->generate_kepler_target(sim_param), target_list, 1:num_sys )
   return KeplerPhysicalCatalog(target_list)
end

function observe_kepler_targets_sky_avg(input::KeplerPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  return observe_kepler_targets(calc_target_obs, input, sim_param)
end

function observe_kepler_targets_single_obs(input::KeplerPhysicalCatalog, sim_param::SimParam )
  calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  return observe_kepler_targets(calc_target_obs, input, sim_param)
end

function observe_kepler_targets(calc_target_obs::Function, input::KeplerPhysicalCatalog, sim_param::SimParam )
  #calc_target_obs = get_function(sim_param,"calc_target_obs_sky_ave")
  #calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
  output = KeplerObsCatalog()
  if haskey(sim_param,"mem_kepler_target_obs")
     output.target = get(sim_param,"mem_kepler_target_obs",Array{KeplerTargetObs}(0) )
  end
  num_targets_sim_pass_one = get_int(sim_param,"num_targets_sim_pass_one")
  if length(output.target) < num_targets_sim_pass_one
     output.target = Array{KeplerTargetObs}(undef,num_targets_sim_pass_one)
  end
  #output.target = Array{KeplerTargetObs}(undef,length(input.target) )  # Replaced to reduce memory allocation
  map!(x::KeplerTarget->calc_target_obs(x,sim_param)::KeplerTargetObs, output.target, input.target)
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
function generate_obs_targets(cat_phys::KeplerPhysicalCatalog, sim_param::SimParam )
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
function simulated_read_kepler_observations(sim_param::SimParam )
   println("# WARNING: Using simulated_read_kepler_observations.")
   # if haskey(sim_param,"stellar_catalog")
   #    star_tab_func = get_function(sim_param, "star_table_setup")
   #    star_tab_func(sim_param)
   # end
   num_sys = get_int(sim_param,"num_kepler_targets")
   generate_kepler_target = get_function(sim_param,"generate_kepler_target")
   target_list = Array{KeplerTarget}(undef,num_sys)
   map!(x->generate_kepler_target(sim_param), target_list, 1:num_sys )

   cat_phys_cut = generate_obs_targets(KeplerPhysicalCatalog(target_list), sim_param)
   calc_target_obs = get_function(sim_param,"calc_target_obs_single_obs")
   output = KeplerObsCatalog()
   output.target = map(x::KeplerTarget->calc_target_obs(x,sim_param)::KeplerTargetObs, cat_phys_cut.target)
   return output
end

function read_koi_catalog(sim_param::SimParam, force_reread::Bool = false)
    filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)),"..", "data", convert(String,get(sim_param,"koi_catalog","q1_q17_dr25_koi.csv")) ) )
    return read_koi_catalog(filename, force_reread)
end

function read_koi_catalog(filename::String, force_reread::Bool = false)
    local df, usable

    if occursin(r".jld2$",filename) && !force_reread
        try
            data = load(filename)
            df = data["koi_catalog"]
            usable = data["koi_catalog_usable"]
            Core.typeassert(df,DataFrame)
            Core.typeassert(usable,Array{Int64,1})
        catch
            error(string("# Failed to read koi catalog >",filename,"< in jld2 format."))
        end
    else
       try
            df = CSV.read(filename,comment="#")

            # Choose which KOIs to keep
            #is_cand = (csv_data[!,:,koi_disposition_idx] .== "CONFIRMED") | (csv_data[!,:,koi_disposition_idx] .== "CANDIDATE")
            is_cand = df[!,:koi_pdisposition] .== "CANDIDATE"
            has_radius = .!ismissing.(df[!,:koi_ror])
            has_period = .!(ismissing.(df[!,:koi_period]) .| ismissing.(df[!,:koi_period_err1]) .| ismissing.(df[!,:koi_period_err2]))

            is_usable = .&(is_cand, has_radius, has_period)
            usable = findall(is_usable)
           #  symbols_to_keep = [:kepid, :kepoi_name, :koi_pdisposition, :koi_score, :koi_ror, :koi_period, :koi_period_err1, :koi_period_err2, :koi_time0bk, :koi_time0bk_err1, :koi_time0bk_err2, :koi_depth, :koi_depth_err1, :koi_depth_err2, :koi_duration, :koi_duration_err1, :koi_duration_err2]
           # df = df[usable, symbols_to_keep]
           # tmp_df = DataFrame()
           # for col in names(df)
           #     tmp_df[col] = collect(skipmissing(df[col]))
           # end
           # df = tmp_df
           # usable = collect(1:length(df[!,:kepid]))
        catch
            error(string("# Failed to read koi catalog >",filename,"< in ascii format."))
        end
    end
    return df, usable
end

# df_star is assumed to have fields kepid, mass and radius for all targets in the survey
function setup_actual_planet_candidate_catalog(df_star::DataFrame, df_koi::DataFrame, usable_koi::Array{Int64}, sim_param::SimParam)
    local target_obs, num_pl
    df_koi = df_koi[usable_koi,:]

    if haskey(sim_param, "koi_subset_csv")
        koi_subset = fill(false, length(df_koi[!,:kepid]))

        subset_df = readtable(convert(String,get(sim_param,"koi_subset_csv", "christiansen_kov.csv")), header=true, separator=' ')

        for n in 1:length(subset_df[!,1])
            subset_colnum = 1
            subset_entry = findall(x->x==subset_df[n,1], df_koi[names(subset_df)[1]])
            # println("Initial cut: ", subset_entry)
            while (length(subset_entry) > 1) & (subset_colnum < length(names(subset_df)))
                subset_colnum += 1

                subsubset = findall(x->round(x*10.)==round(subset_df[n,subset_colnum]*10.), df_koi[subset_entry,names(subset_df)[subset_colnum]])
	        # println("Extra cut: ", subset_df[n,subset_colnum], " / ", df_koi[subset_entry,col_idx], " = ", subsubset)
	        subset_entry = subset_entry[subsubset]
            end
            if length(subset_entry) > 1
	        cand_sub = findall(x->x == "CANDIDATE",df_koi[subset_entry,:koi_pdisposition])
	        subset_entry = subset_entry[cand_sub]
	        if length(subset_entry) > 1
                    println("# Multiple planets found in final cut: ", subset_df[n,1])
                end
            end
            if length(subset_entry) < 1
                println("# No planets found in final cut: ", subset_df[n,:])
            end
            koi_subset[subset_entry] = true
        end
        df_koi = df_koi[findall(koi_subset),:]
        tot_plan = count(x->x, koi_subset)
    end

    output = KeplerObsCatalog()
    sort!(df_star, (:kepid))
    df_obs = join(df_star, df_koi, on = :kepid)
    #df_obs = sort!(df_obs, cols=(:kepid))
    df_obs = sort!(df_obs, (:kepid))

    if haskey(sim_param, "koi_subset_csv")
        tot_plan -= length(df_obs[!,:kepoi_name])
        println("# Number of planet candidates in subset file with no matching star in table: ", tot_plan)
    end

    plid = 0
    for i in 1:length(df_obs[!,:kepoi_name])
        if plid == 0
            plid = 1
            while i+plid < length(df_obs[!,:kepoi_name]) && df_obs[i+plid,:kepid] == df_obs[i,:kepid]
                plid += 1
            end
            num_pl = plid
            target_obs = KeplerTargetObs(num_pl)
	        #target_obs.star = ExoplanetsSysSim.StarObs(df_obs[i,:radius],df_obs[i,:mass],findfirst(df_star[!,:kepid], df_obs[i,:kepid]))
            star_idx = searchsortedfirst(df_star[!,:kepid],df_obs[i,:kepid])
            if star_idx > length(df_star[!,:kepid])
                @warn "# Couldn't find kepid " * df_star[i,:kepid] * " in df_obs."
                star_idx = rand(1:length(df_star[!,:kepid]))
            end
            target_obs.star = ExoplanetsSysSim.StarObs(df_obs[i,:radius],df_obs[i,:mass],star_idx)

        end

        target_obs.obs[plid] = ExoplanetsSysSim.TransitPlanetObs(df_obs[i,:koi_period],df_obs[i,:koi_time0bk],df_obs[i,:koi_depth]/1.0e6,df_obs[i,:koi_duration])
        target_obs.sigma[plid] = ExoplanetsSysSim.TransitPlanetObs((abs(df_obs[i,:koi_period_err1])+abs(df_obs[i,:koi_period_err2]))/2,(abs(df_obs[i,:koi_time0bk_err1])+abs(df_obs[i,:koi_time0bk_err2]))/2,(abs(df_obs[i,:koi_depth_err1]/1.0e6)+abs(df_obs[i,:koi_depth_err2]/1.0e6))/2,(abs(df_obs[i,:koi_duration_err1])+abs(df_obs[i,:koi_duration_err2]))/2)
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
function calc_snr_list(cat::KeplerPhysicalCatalog, sim_param::SimParam)
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

function calc_prob_detect_list(cat::KeplerPhysicalCatalog, sim_param::SimParam)
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
  cat_phys = generate_kepler_physical_catalog(sim_param)::KeplerPhysicalCatalog
  id = findfirst( x->num_planets(x)>=1 , cat_phys.target)   # fast forward to first target that has some planets
  @assert(length(id)>=1)
  semimajor_axis(cat_phys.target[id].sys[1],1)
  pdetlist = calc_prob_detect_list(cat_phys,sim_param)
  calc_target_obs_single_obs(cat_phys.target[id],sim_param)
  calc_target_obs_sky_ave(cat_phys.target[id],sim_param)
  @assert( length(cat_phys.target[id].sys[1].planet)  == num_planets(cat_phys.target[id]) )
  cat_obs = simulated_read_kepler_observations(sim_param)
  return (cat_phys, cat_obs)
end
