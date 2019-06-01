## ExoplanetsSysSim.jl
## (c) 2015 Eric B. Ford

module SimulationParameters

import Compat: @compat #, readstring
import Pkg, LibGit2
using ExoplanetsSysSim

export SimParam, add_param_fixed, add_param_active, update_param, set_active, set_inactive, is_active
export get_any, get_real, get_int, get_bool, get_function

export make_vector_of_active_param_keys, make_vector_of_sim_param, get_range_for_sim_param, update_sim_param_from_vector!
export setup_sim_param_demo, test_sim_param_constructors
#export preallocate_memory!

const global version_id_str = "1.0.0"
const global version_id_pair = ("version",version_id_str)
const global julia_version_pair = ("version_julia",string(VERSION))

mutable struct SimParam
  param::Dict{String,Any}
  active::Dict{String,Bool}
end
copy(p::SimParam) = SimParam(copy(p.param),copy(p.active))

"""
    SimParam(p::Dict{String,Any})
Creates a SimParam object from the dictionary p, with all parameter defaulting to inactive.
"""
function SimParam(p::Dict{String,Any})   # By default all parameters are set as inactive (i.e., not allowed to be optimized)
  a = Dict{String,Bool}()
  for k in keys(p)
    a[k] = false
  end
  return SimParam(p,a)
end

function SimParam()
  d = Dict{String,Any}([version_id_pair,julia_version_pair])
  return SimParam(d)
end

"Update SimParam() to define current state at run time, and not precompile time"
function __init__()
"""
    SimParam()
Creates a nearly empty SimParam object, with just the version id and potentially other information about the code, system, runtime, etc.
"""
function SimParam()
  d = Dict{String,Any}([version_id_pair,julia_version_pair,("hostname",gethostname()), ("time",time())])
  try
    d["ExoplanetsSysSim version"] = Pkg.installed()["ExoplanetsSysSim"]
    #d["ExoplanetsSysSim directory"] = dirname(pathof(ExoplanetsSysSim))
    #d["ExoplanetsSysSim branch"] = LibGit2.headname(LibGit2.GitRepo(dirname(pathof(ExoplanetsSysSim))))
    #d["ExoplanetsSysSim head_oid"] = LibGit2.head_oid(LibGit2.GitRepo(dirname(pathof(ExoplanetsSysSim))))
  catch
    @warn("# Couldn't add full information about version of SysSim to SimParam.")
  end
  SimParam(d)
end
end

"""
    add_param_fixed(sim::SimParam, key::String,val::Any)
Adds (or overwrites) key with value val to the SimParam object, sim, and sets the parameter set to inactive.
"""
function add_param_fixed(sim::SimParam, key::String,val::Any)
  sim.param[key] = val
  sim.active[key] = false
end

"""
### add_param_active(sim::SimParam, key::String,val::Any)
Adds (or overwrites) key with value val to the SimParam object, sim, and sets the parameter set to active.
"""
function add_param_active(sim::SimParam, key::String,val::Any)
  sim.param[key] = val
  sim.active[key] = true
end

"""
### update_param(sim::SimParam, key::String,val::Any)
Overwrites key with value val to the SimParam object, sim
"""
function update_param(sim::SimParam, key::String,val::Any)
  @assert haskey(sim.param,key)
  sim.param[key] = val
end

"""
### set_active(sim::SimParam, key::String)
Sets the key parameter to be active in sim.
"""
function set_active(sim::SimParam,key::String)
  @assert haskey(sim.param,key)
  sim.active[key] = true
end

"""
### set_active(sim::SimParam, keys::Vector{String})
Sets each of the key parameters to be active in sim.
"""
function set_active(sim::SimParam,keys::Vector{String})
  for k in keys
    set_active(sim,k)
  end
end

"""
### set_inactive(sim::SimParam, key::String)
Sets the key parameter to be inactive in sim.
"""
function set_inactive(sim::SimParam,key::String)
  @assert haskey(sim.param,key)
  sim.active[key] = false
end

"""
### set_inactive(sim::SimParam, keys::Vector{String})
Sets each of the key parameters to be inactive in sim.
"""
function set_inactive(sim::SimParam,keys::Vector{String})
  for k in keys
    set_inactive(sim,k)
  end
end

function is_active(sim::SimParam,key::String)
  @assert haskey(sim.active,key)
  sim.active[key]
end

import Base.get
function get(sim::SimParam, key::String, default_val::T) where T
  val::T = get(sim.param,key,default_val)::T
  return val
end

function get_any(sim::SimParam, key::String, default_val::Any)
  val = get(sim.param,key,default_val)
  return val
end

function get_real(sim::SimParam, key::String)
  val::Float64 = get(sim.param,key,convert(Float64,NaN) )::Float64
  @assert(val!=convert(Float64,NaN))
  return val
end

function get_int(sim::SimParam, key::String)
  val::Int64 = get(sim.param,key,zero(Int64))
  #@assert(val!=nan(zero(Int64)))
  #@assert(val!=oftype(x,NaN))
  return val
end

function get_bool(sim::SimParam, key::String)
  val::bool = get(sim.param,key,false)
  return val
end

function noop()
end


function get_function(sim::SimParam, key::String)
  val::Function = Base.get(sim.param,key,noop)::Function
  #val = Base.get(sim.param,key,null)
  @assert((val!=nothing) && (val!=noop))
  return val::Function
end

import Base.haskey
haskey(sim::SimParam, key::String) = return haskey(sim.param,key)

function make_vector_of_active_param_keys(sim::SimParam)
  sortedkeys = sort(collect(keys(sim.param)))
  sortedkeys_active = sortedkeys[map(k->get(sim.active,k,false),sortedkeys)]
  return sortedkeys_active
end


function make_vector_of_sim_param(sim::SimParam)
  param_vector = Float64[]      # QUERY: Currently we make separate vectors of Ints and Floats.  Does this make sense?
  sss = sort(collect(keys(sim.param)))
  for k in 1:length(sss)
    if(sim.active[sss[k]]==false)
      continue
    end
    if(length(sim.param[sss[k]])==1)
      if isa( sim.param[sss[k]], Real )
        push!(param_vector,sim.param[sss[k]])
      elseif(eltype( sim.param[sss[k]]) <: Real)
        append!(param_vector,vec(sim.param[sss[k]]))
      end
    elseif(eltype( sim.param[sss[k]]) <: Real)
      append!(param_vector, vec(sim.param[sss[k]]))
    else
      if eltype( sss[k]) <: Real
        append!(param_vector,sim.param[sss[k]])
      end
    end
  end
  return param_vector
end

function get_range_for_sim_param(key::String, sim::SimParam)
  sorted_keys = sort(collect(keys(sim.param)))
  i = 1
  for k in 1:length(sorted_keys)
    if(sim.active[sorted_keys[k]]==false)
      continue
    end
    param_len = length(sim.param[sorted_keys[k]])
    if sorted_keys[k]==key
       return i:(i+param_len)
    else
       #println("Didn't match >",sorted_keys[k],"< and >",key,"<.")
       i += param_len
    end
  end
  println("# ERROR: Never found range for param: ",key)
  return 0:0
end

function update_sim_param_from_vector!(param::Vector{Float64}, sim::SimParam)
  #println("# Input vector: ",param)
  sorted_keys = sort(collect(keys(sim.param)))
  i = 1
  for k in 1:length(sorted_keys)
    if(sim.active[sorted_keys[k]]==false)
      continue
    end
    param_len = length(sim.param[sorted_keys[k]])
    if param_len==1
      if isa( sim.param[sorted_keys[k]], Real )
        # println("# Replacing >",sorted_keys[k],"< with >",param[i],"<")
        sim.param[sorted_keys[k]] = param[i]
        i = i+1
      elseif eltype( sim.param[sorted_keys[k]]) <: Real
        # println("# Replacing >",sim.param[sorted_keys[k]],"< with >",reshape(param[i:i+param_len-1], size(sim.param[sorted_keys[k]]),"<")
        sim.param[sorted_keys[k]] = reshape(param[i:i+param_len-1], size(sim.param[sorted_keys[k]]))
        i = i+1
      end
    elseif param_len>1
      if eltype( sim.param[sorted_keys[k]]) <: Real
        # println("# Replacing >",sim.param[sorted_keys[k]],"< with >",reshape(param[i:i+param_len-1], size(sim.param[sorted_keys[k]])),"<")
        sim.param[sorted_keys[k]] = reshape(param[i:i+param_len-1], size(sim.param[sorted_keys[k]]))
        i = i+param_len
      end
    else
      println("# Don't know what to do with empty simulation parameter: ",sorted_keys[k])
    end
  end
  return sim
end

function preallocate_memory!(sim_param::SimParam)
  num_kepler_targets = get_int(sim_param,"num_kepler_targets")
  add_param_fixed(sim_param,"mem_kepler_target_obs", Array{KeplerTargetObs}(num_kepler_targets) )
end

function setup_sim_param_demo(args::Vector{String} = Array{String}(undef,0) )   # allow this to take a list of parameter (e.g., from command line)
  sim_param = SimParam()
  add_param_fixed(sim_param,"max_tranets_in_sys",7)
  add_param_fixed(sim_param,"num_targets_sim_pass_one",190000)                      # Note this is used for the number of stars in the simulations, not necessarily related to number of Kepler targets
  add_param_fixed(sim_param,"num_kepler_targets",190000)                           # Note this is used for the number of Kepler targets for the observational catalog
  add_param_fixed(sim_param,"generate_star",ExoplanetsSysSim.generate_star_dumb)
  #add_param_fixed(sim_param,"generate_planetary_system", ExoplanetsSysSim.generate_planetary_system_simple)
  add_param_fixed(sim_param,"generate_planetary_system", ExoplanetsSysSim.generate_planetary_system_uncorrelated_incl)

  # add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_simple)
  add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_from_table)
  add_param_fixed(sim_param,"star_table_setup",StellarTable.setup_star_table)
  add_param_fixed(sim_param,"stellar_catalog","q1q17_dr25_gaia_fgk.jld2")
  add_param_fixed(sim_param,"generate_num_planets",ExoplanetsSysSim.generate_num_planets_poisson)
  add_param_active(sim_param,"log_eta_pl",log(2.0))
  add_param_fixed(sim_param,"generate_planet_mass_from_radius",ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw)
  add_param_fixed(sim_param,"vetting_efficiency",ExoplanetsSysSim.vetting_efficiency_none)
  add_param_fixed(sim_param,"mr_power_index",2.0)
  add_param_fixed(sim_param,"mr_const",1.0)
  #add_param_fixed(sim_param,"generate_period_and_sizes",ExoplanetsSysSim.generate_period_and_sizes_log_normal)
  #add_param_active(sim_param,"mean_log_planet_radius",log(2.0*earth_radius))
  #add_param_active(sim_param,"sigma_log_planet_radius",log(2.0))
  #add_param_active(sim_param,"mean_log_planet_period",log(5.0))
  #add_param_active(sim_param,"sigma_log_planet_period",log(2.0))
  add_param_fixed(sim_param,"generate_period_and_sizes", ExoplanetsSysSim.generate_period_and_sizes_power_law)
  add_param_active(sim_param,"power_law_P",0.3)
  add_param_active(sim_param,"power_law_r",-2.44)
  add_param_fixed(sim_param,"min_period",1.0)
  add_param_fixed(sim_param,"max_period",100.0)
  add_param_fixed(sim_param,"min_radius",0.5*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"max_radius",10.0*ExoplanetsSysSim.earth_radius)
  add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
  add_param_fixed(sim_param,"sigma_hk",0.03)
  add_param_fixed(sim_param,"sigma_incl",2.0)   # degrees
  add_param_fixed(sim_param,"calc_target_obs_sky_ave",ExoplanetsSysSim.calc_target_obs_sky_ave)
  add_param_fixed(sim_param,"calc_target_obs_single_obs",ExoplanetsSysSim.calc_target_obs_single_obs)
  add_param_fixed(sim_param,"read_target_obs",ExoplanetsSysSim.simulated_read_kepler_observations)
  add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_fixed_noise)
  # add_param_fixed(sim_param,"transit_noise_model",transit_noise_model_diagonal)
  # add_param_fixed(sim_param,"rng_seed",1234)   # If you want to be able to reproduce simulations


  # Do other initialization tasks belong here or elsewhere?
  # TODO OPT:  Try to preallocate memory for each target to see if this makes a performance difference
  # preallocate_memory!(sim_param)

  return sim_param
end


function test_sim_param_constructors()
  oldval = log(2.0)
  sim_param = SimParam( Dict([ ("version",version_id_str), ("num_kepler_targets",190000), ("log_eta_pl",oldval), ("max_tranets_in_sys",7)] ) )
  get(sim_param,"version","")
  set_active(sim_param,"log_eta_pl")
  sp_vec = make_vector_of_sim_param(sim_param)
  sp_vec .+= 0.1
  update_sim_param_from_vector!(sp_vec,sim_param)
  newval = get_real(sim_param,"log_eta_pl")
  isapprox(oldval+0.1,newval,atol=0.001)
end

end

#test_sim_param_constructors()
