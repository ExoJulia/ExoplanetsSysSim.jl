## ExoplanetsSysSim/src/stellar_table.jl
## (c) 2015 Eric B. Ford

module StellarTable
using ExoplanetsSysSim
#using DataArrays
using DataFrames
using CSV
using JLD
using JLD2
using FileIO

#if VERSION >= v"0.5-"
#  import Compat: UTF8String, ASCIIString
#end

export setup_star_table, star_table, num_usable_in_star_table, set_star_table, star_table_has_key

df = DataFrame()

function setup(sim_param::SimParam; force_reread::Bool = false)
  global df
  wf = WindowFunction.setup_window_function(sim_param)
  WindowFunction.setup_OSD_interp(sim_param) #read in osd files so they can be interpolated
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)), "..","data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr24_stellar.csv")) ) )
  df = setup(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  add_param_fixed(sim_param,"num_kepler_targets",num_usable_in_star_table())
  return df
end

function setup(filename::String; force_reread::Bool = false)
  global df
  if occursin(r".jld2$",filename)
  #if occursin(r".jld$",filename)
  try
    data = load(filename)
    df = data["stellar_catalog"]
    Core.typeassert(df,DataFrame)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld2 format."))
  end
  else
  try
    #df = readtable(filename)
    #df = CSV.read(filename,nullable=true)
    df = CSV.read(filename, allowmissing=:all)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  # Now we read in all CDPP's, so can interpolate to transit duration
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :cdppslplong, :cdppslpshrt, :dataspan, :dutycycle, :limbdark_coeff1, :limbdark_coeff2, :limbdark_coeff3, :limbdark_coeff4 ]

  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  is_usable = [ !any(ismissing.([ df[i,j] for j in 1:size(df,2) ])) for i in 1:size(df,1) ]
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
  end
    df[:wf_id] = map(x->ExoplanetsSysSim.WindowFunction.get_window_function_id(x,use_default_for_unknown=false),df[:kepid])
    obs_5q = df[:wf_id].!=-1
    #df = df[obs_5q,keys(df.colindex)]
    df = df[obs_5q,names(df)]
    StellarTable.set_star_table(df)
  return df
end

setup_star_table(sim_param::SimParam; force_reread::Bool = false) = setup(sim_param, force_reread=force_reread)
setup_star_table(filename::String) = setup(filename)

function num_usable_in_star_table()
  global df
  return size(df,1)
end

function star_table(i::Integer, sym::Symbol)
  global df
  return df[i,sym]
end

function star_table(i::Integer)
  global df
  return df[i,:]
end

function star_table(i::Integer, sym::Vector{Symbol})
  global df
  return df[i,sym]
end

function star_table(i::Vector{Integer}, sym::Symbol)
  global df
  return df[i,sym]
end

function star_table(i::Vector{Integer}, sym::Vector{Symbol})
  global df
  return df[i,sym]
end

function set_star_table(df2::DataFrame)
  global df
  df = df2
end

function star_table_has_key(s::Symbol)
  global df
  haskey(df,s)
end

end # module StellarTable

# using ExoplanetsSysSim.StellarTable

#=
function generate_star_from_table(sim_param::SimParam, id::Integer)  # WARNING:  To be renamed once there's a working/tested version that uses a stellar catalog with GAIA data
  mu_r = StellarTable.star_table(id,:radius)
  sig_r1 = StellarTable.star_table(id,:radius_err1)
  sig_r2 = StellarTable.star_table(id,:radius_err2)
  z = randn()
  r = mu_r + (z>0) ?  z*sig_r1 : z*sig_r2
  m = rand(Normal(r,0.1))::Float64
  while m<0.0
    m = rand(Normal(r,0.1))::Float64
  end
  f = rand(Normal(1.0,0.1))::Float64
  while f<0.0
    f = 1.0+0.1*randn()
  end
  ld = LimbDarkeningParam4thOrder(StellarTable.star_table(id,:limbdark_coeff1), StellarTable.star_table(id,:limbdark_coeff2), StellarTable.star_table(id,:limbdark_coeff3), StellarTable.star_table(id,:limbdark_coeff4) )
  return SingleStar(r,m,f,id,ld)
  #return SingleStar(r,m,f,id)
end
=#

function generate_star_from_table(sim_param::SimParam)
  id = rand(1:StellarTable.num_usable_in_star_table())
  generate_star_from_table(sim_param, id)
end
