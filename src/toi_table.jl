## ExoplanetsSysSim/src/toi_table.jl
## (c) 2015 Eric B. Ford

# Note this file is currently not used by SysSim.
# This functionality is now in kepler_catalog.jl

module KoiTable
using ExoplanetsSysSim
#using DataArrays
using DataFrames
using CSV

export setup_toi_table, toi_table, num_toi_for_ticid

df = DataFrame()
usable = Array{Int64}(0)
        
default_toi_symbols_to_keep = [ :toi_name, :toi_vet_stat, :toi_pdisposition, :toi_period, :toi_time0bk, :toi_duration, :toi_ingress, :toi_depth, :toi_ror, :toi_prad, :toi_srad, :toi_smass, :toi_steff, :toi_slogg, :toi_smet ]

function setup(sim_param::SimParam; force_reread::Bool = false, symbols_to_keep::Vector{Symbol} = default_toi_symbols_to_keep )
  global df, usable
  if haskey(sim_param,"read_toi_catalog") && !force_reread
     return df
  end
  toi_catalog = joinpath(dirname(pathof(ExoplanetsSysSim)),"..", "data", get(sim_param,"toi_catalog","tois.csv") )
  add_param_fixed(sim_param,"read_toi_catalog",true)
  try 
    #df = readtable(toi_catalog)
    #df = CSV.read(toi_catalog,nullable=true)
    df = CSV.read(toi_catalog,allowmissing=:all)
  catch
    error(string("# Failed to read toi catalog >",toi_catalog,"<."))
  end

  has_planet = ! (isna(df[:toi_period]) | isna(df[:toi_time0bk]) | isna(df[:toi_duration]) | isna(:toi_depth) )
  has_star = ! ( isna(:toi_srad) )
  is_usable = has_planet & has_star

  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = find(is_usable)
  df = df[usable, symbols_to_keep]
end

setup_toi_table(sim_param::SimParam) = setup(sim_param::SimParam)

function ticids_w_tois()
  unique(df[:,:ticid])
end

function df_for_ticid(ticid::Integer)
  df[df[:ticid].==ticid,:]
end

function num_toi(ticid::Integer)
  sum(df[:ticid].==ticid)
end

function toi_by_ticid(ticid::Integer, plid::Integer, sym::Symbol)
  ticid_idx = df[:ticid].==ticid
  per_perm = sortperm(df[ticid_idx,:toi_period])
  @assert( 1<= plid <= length(per_perm) )
  df[ticid_idx,sym][per_perm[plid]]
end


function num_usable()
  global usable
  length(usable)
end

num_usable_in_toi_table() = num_usable()

function idx(i::Integer)
  global usable
  @assert( 1<=i<=length(usable) )
  usable[i]
end


function toi_table(i::Integer, sym::Symbol)
  global df, usable
  @assert( 1<=i<=length(usable) )
  return df[i,sym]
  #return df[usable[i],sym]
end

function toi_table(i::Integer)
  global data
  return df[i,:]
  #return df[usable[i],:]
end

function toi_table(i::Integer, sym::Vector{Symbol})
  global df, usable
  @assert( 1<=i<=length(usable) )
  return df[i,sym]
  #return df[usable[i],sym]
end

function toi_table(i::Vector{Integer}, sym::Symbol)
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

function toi_table(i::Vector{Integer}, sym::Vector{Symbol})
  global df, usable
  return df[i,sym]
  #return df[usable[i],sym]
end

end # module KoiTable

using ExoplanetsSysSim.KoiTable

