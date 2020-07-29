## ExoplanetsSysSim/src/koi_table.jl
## (c) 2015 Eric B. Ford

# Note this file is currently not used by SysSim.
# This functionality is now in kepler_catalog.jl

module KoiTable
using ExoplanetsSysSim
#using DataArrays
using DataFrames
using CSV

export setup_koi_table, koi_table, num_koi_for_kepid

df = DataFrame()
usable = Array{Int64}(0)
        
default_koi_symbols_to_keep = [ :kepoi_name, :koi_vet_stat, :koi_pdisposition, :koi_period, :koi_time0bk, :koi_duration, :koi_ingress, :koi_depth, :koi_ror, :koi_prad, :koi_srad, :koi_smass, :koi_steff, :koi_slogg, :koi_smet ]

function setup(sim_param::SimParam; force_reread::Bool = false, symbols_to_keep::Vector{Symbol} = default_koi_symbols_to_keep )
  global df, usable
  if haskey(sim_param,"read_koi_catalog") && !force_reread
     return df
  end
  planetary_catalog = joinpath(dirname(pathof(ExoplanetsSysSim)),"..", "data", get(sim_param,"planetary_catalog","q1_q17_dr25_koi.csv") )
  add_param_fixed(sim_param,"read_koi_catalog",true)
  try 
    #df = readtable(planetary_catalog)
    #df = CSV.read(planetary_catalog,nullable=true)
    df = CSV.read(planetary_catalog,allowmissing=:all)
  catch
    error(string("# Failed to read koi catalog >",planetary_catalog,"<."))
  end

  has_planet = ! (isna(df[:koi_period]) | isna(df[:koi_time0bk]) | isna(df[:koi_duration]) | isna(:koi_depth) )
  has_star = ! ( isna(:koi_srad) )
  is_usable = has_planet & has_star

  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
end

setup_koi_table(sim_param::SimParam) = setup(sim_param::SimParam)

function kepids_w_kois()
  unique(df[:,:kepid])
end

function df_for_kepid(kepid::Integer)
  df[df[:kepid].==kepid,:]
end

function num_koi(kepid::Integer)
  sum(df[:kepid].==kepid)
end

function koi_by_kepid(kepid::Integer, plid::Integer, sym::Symbol)
  kepid_idx = df[:kepid].==kepid
  per_perm = sortperm(df[kepid_idx,:koi_period])
  @assert( 1<= plid <= length(per_perm) )
  df[kepid_idx,sym][per_perm[plid]]
end


function num_usable()
  global usable
  length(usable)
end

num_usable_in_koi_table() = num_usable()

function idx(i::Integer)
  global usable
  @assert( 1<=i<=length(usable) )
  usable[i]
end


function koi_table(i::Integer, sym::Symbol)
  global df, usable
  @assert( 1<=i<=length(usable) )
  return df[i,sym]
  #return df[usable[i],sym]
end

function koi_table(i::Integer)
  global data
  return df[i,:]
  #return df[usable[i],:]
end

function koi_table(i::Integer, sym::Vector{Symbol})
  global df, usable
  @assert( 1<=i<=length(usable) )
  return df[i,sym]
  #return df[usable[i],sym]
end

function koi_table(i::Vector{Integer}, sym::Symbol)
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

function koi_table(i::Vector{Integer}, sym::Vector{Symbol})
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

end # module KoiTable

using ExoplanetsSysSim.KoiTable

