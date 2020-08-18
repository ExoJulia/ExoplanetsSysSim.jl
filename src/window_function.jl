## ExoplanetsSysSim/src/window_function.jl
## (c) 2018 Darin Ragozzine

# Gather and prepare the window function data
module WindowFunction

import ..cdpp_durations
export setup_window_function, get_window_function_data, get_window_function_id, eval_window_function

#using DataArrays
using DataFrames
#using CSV
#using JLD2
using FileIO
using ExoplanetsSysSim
using ExoplanetsSysSim.SimulationParameters


# Object to hold window function data
struct window_function_data
  window_func_array::Array{Float64,3}      #  Value of window function (window_function_id, duration_id, period_id).  Maybe rename to wf_value or data?
  wf_durations_in_hrs::Array{Float64,1}    # Boundaries for duration bins in window_func_array.  Maybe rename to durations?
  wf_periods_in_days::Array{Float64,1}     # Boundaries for periods bins in window_func_array.  Maybe rename to periods?
  sorted_quarter_strings::Array{Int64,1}   # TODO OPT: Is there a reason to keep this?  Maybe rename to quarter_strings?
  allsortedkepids::Array{Int64,1}          # value is Kepler ID.  Index is same as index to window_function_id_arr
  window_function_id_arr::Array{Int64,1}   # value is index for window_func_array.  Index is same as index to allsortedkepids
  default_wf_id::Int64                     # Index corresponding to the default window function
  already_warned::Array{Bool,1}            # Wether we've already thrown a warning about this kepid
end



function window_function_data()
  window_function_data( Array{Float64,3}(undef,0,0,0), Array{Float64,1}(undef,0),Array{Float64,1}(undef,0), Array{Int64,1}(undef,0),Array{Int64,1}(undef,0),Array{Int64,1}(undef,0), 0, falses(0) )
end


win_func_data = window_function_data()

function setup(sim_param::SimParam; force_reread::Bool = false)
  global win_func_data
  if haskey(sim_param,"read_window_function") && !force_reread
     return win_func_data
  end
  window_function_filename = convert(String,joinpath(dirname(pathof(ExoplanetsSysSim)),
      "..", "data", convert(String,get(sim_param,"window_function","DR25topwinfuncs.jld2")) ) )
  setup(window_function_filename)
  add_param_fixed(sim_param,"read_window_function",true)
  @assert( size(win_func_data.window_func_array,2) == length(win_func_data.wf_durations_in_hrs) )
  @assert( size(win_func_data.window_func_array,3) == length(win_func_data.wf_periods_in_days) )
  @assert( size(win_func_data.window_func_array,1) >= maximum(win_func_data.window_function_id_arr) )
  @assert( size(win_func_data.window_func_array,1) >= win_func_data.default_wf_id )
  return win_func_data
end



function setup(filename::String)
# Reads in the window function data collected from the Kepler Completeness Products
# see Darin Ragozzine's get/cleanDR25winfuncs.jl

  if occursin(r".jld2$",filename)
    try
      wfdata = load(filename)
      window_func_array = wfdata["window_func_array"]
      wf_durations_in_hrs = wfdata["wf_durations_in_hrs"]  # TODO OPT DETAIL: Should we convert units to days here?
      wf_periods_in_days = wfdata["wf_periods_in_days"]
      sorted_quarter_strings = wfdata["sorted_quarter_strings"]
      allsortedkepids = wfdata["allsortedkepids"]
      window_function_id_arr = wfdata["window_function_id_arr"]
      already_warned = falses(length(allsortedkepids))
      global win_func_data = window_function_data(window_func_array, wf_durations_in_hrs, wf_periods_in_days, sorted_quarter_strings,
                                                allsortedkepids, window_function_id_arr, maximum(window_function_id_arr), already_warned )

    catch
      error(string("# Failed to read window function data > ", filename," < in jld2 format."))
    end
  end

  return win_func_data
end

setup_window_function(sim_param::SimParam; force_reread::Bool = false) = setup(sim_param, force_reread=force_reread)
setup_window_function(filename::String; force_reread::Bool = false) = setup(filename, force_reread=force_reread)

function get_window_function_data()::window_function_data
   #global win_func_data
   return win_func_data
end

function get_window_function_id(kepid::Int64; use_default_for_unknown::Bool = true)::Int64
  # takes the quarter string from the stellar catalog and determines the window function id
  # from DR25topwinfuncs.jld2 made by Darin Ragozzine's cleanDR25winfuncs.jl script.
  no_win_func_available::Int64 = -1        # hardcoding this in, should match convention in window function input file

  idx = searchsortedfirst(win_func_data.allsortedkepids,kepid) # all Kepler kepids are in allsortedkepids
  wf_id = win_func_data.window_function_id_arr[idx]


  if wf_id == no_win_func_available && use_default_for_unknown
    # if a target is observed for less than 4 quarters, then it won't have a corresponding
    # window function in this list, so throw a warning and use the last window_function_id
    # which corresponds to an "averaged" window function
    if !win_func_data.already_warned[idx]
       win_func_data.already_warned[idx] = true
       if sum(win_func_data.already_warned) < 20
          @warn "Window function data is not avaialble for kepid $kepid, using default."
       end
    end
    wf_id = win_func_data.default_wf_id
  end
  # TODO SCI DETAIL IMPORTANT? This does not include TPS timeouts or MESthresholds (see DR25 Completeness Products)

  return wf_id
end


function calc_period_idx(P::Float64)::Int64
  @assert(P>zero(P))
  idx = searchsortedlast(win_func_data.wf_periods_in_days,P)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_periods_in_days)
     if P-win_func_data.wf_periods_in_days[idx]>win_func_data.wf_periods_in_days[idx+1]-P
        idx += 1
     end
  end
  return idx
end

function calc_duration_idx(D::Float64)::Int64
  # NOTE: Currently assumes we left wf data in hours, so deal with that conversion here
  @assert(D>=zero(D)) ##### Make sure this function is still doing the right thing if D = 0!
  hours_in_day = 24
  idx = searchsortedlast(win_func_data.wf_durations_in_hrs,D*hours_in_day)
  if idx == 0
     return 1
  elseif idx<length(win_func_data.wf_durations_in_hrs)
     if D*hours_in_day-win_func_data.wf_durations_in_hrs[idx]>win_func_data.wf_durations_in_hrs[idx+1]-D*hours_in_day
        idx += 1
     end
  end
  return idx
end


function eval_window_function(wf_idx::Int64=-1; Duration::Float64=0., Period::Float64=0.)::Float64
  D_idx = calc_duration_idx(Duration)
  P_idx = calc_period_idx(Period)
  wf = eval_window_function(wf_idx,D_idx,P_idx)
  # TODO SCI DETAIL: Improve way deal with missing wf values for some durations. Interpolate?
  while wf<=zero(wf) && D_idx<length(win_func_data.wf_durations_in_hrs)
     D_idx += 1
     wf = eval_window_function(wf_idx,D_idx,P_idx)
  end
  return wf
end

function eval_window_function(wf_idx::Int64, D_idx::Int64, P_idx::Int64)::Float64
   global win_func_data
   #@assert(1<=wf_idx<maximum(win_func_data.window_function_id_arr))
   #@assert(1<=P_idx<=length(win_func_data.wf_periods_in_days))
   #@assert(1<=D_idx<=length(win_func_data.wf_durations_in_hrs))
   return win_func_data.window_func_array[wf_idx,D_idx,P_idx]
end

#Object for storing data necessary for OSD_interpolator
struct OSD_data{T1<:Real, T2<:Real}
    allosds::Array{T1,3}
    kepids::Array{Int64,1}
    #periods_length::Int64
    #durations_length::Int64
    grid::Array{Array{T2,1},1}

    function OSD_data(data::AbstractArray{T1,3}, kepids::AbstractArray{Int64,1}, durations::AbstractArray{T2,1}, periods::AbstractArray{T2,1} ) where {T1<:Real, T2<:Real}
        @assert(size(data,1)==length(kepids))
        @assert(size(data,2)==length(durations))
        @assert(size(data,3)==length(periods))
        @assert issorted(kepids)
        @assert issorted(durations)
        @assert issorted(periods)
        new{T1,T2}(data, kepids, [durations,periods])
    end
end

# Only point of this version is to provide a drop in replacement for Keir's code
function OSD_data(allosds::AbstractArray{T1,3}, kepids::AbstractArray{T3,1}, periods_length::Int64, durations_length::Int64, grid::Array{Array{T2,1},1}) where {T1<:Real, T2<:Real, T3<:Real}
    @assert length(grid) == 2
    @assert durations_length == length(grid[1])
    @assert periods_length == length(grid[2])
    if eltype(kepids) != Int64
        kepids = round.(Int64,kepids)
    end
    @assert grid[1][1] == 1.5  # Checking that durations were passed as hours, since that's what Keir's assumed
    return OSD_data(allosds, kepids, grid[1] ./ 24.0 , grid[2])  # Convert durations to days, since that's units in rest of SysSim
end

num_stars(osd::OSD_data) = size(osd.allosds,1)
num_durations(osd::OSD_data) = size(osd.allosds,2)
num_periods(osd::OSD_data) = size(osd.allosds,3)

function setup_OSD(sim_param::SimParam; force_reread::Bool = false)			#reads in 3D table of OSD values and sets up global variables to be used in interpolation
  global OSD_setup
  if haskey(sim_param,"read_OSD_function") && !force_reread
     return OSD_setup
  end
  #OSD_file = load(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", convert(String,get(sim_param,"osd_file","allosds.jld"))))
  #OSD_file = load(joinpath(Pkg.dir(), "ExoplanetsSysSim", "data", convert(String,get(sim_param,"osd_file","allosds.jld"))))
  #OSD_file = load(joinpath(dirname(pathof(ExoplanetsSysSim)),"data",convert(String,get(sim_param,"osd_file","allosds.jld"))))
  #OSD_file = load(joinpath(dirname(pathof(ExoplanetsSysSim)),"..","data",convert(String,get(sim_param,"osd_file","dr25fgk_relaxcut_osds.jld"))))
  OSD_file = load(joinpath(dirname(pathof(ExoplanetsSysSim)),"..","data",convert(String,get(sim_param,"osd_file","dr25fgk_small_osds.jld2"))))
  allosds = OSD_file["allosds"]			#table of OSDs with dimensions: kepids,durations,periods
  periods = OSD_file["periods"][1,:]		#1000 period values corresponding to OSD values in the third dimension of the allosds table
  kepids = OSD_file["kepids"]			#kepids corresponding to OSD values in the first dimension of the allosds table
  OSD_file = 0 # unload OSD file to save memory
  #durations = [1.5,2.,2.5,3.,3.5,4.5,5.,6.,7.5,9.,10.5,12.,12.5,15.] #14 durations corresponding to OSD values in the first dimension of theh allosds table
  periods_length = length(allosds[1,1,:])
  durations_length = length(allosds[1,:,1])
  @assert length(cdpp_durations) == durations_length
  grid = Array{Float64,1}[]			#grid used in OSD_interpolator
  push!(grid, cdpp_durations)
  push!(grid, periods)
  #global compareNoise = Float64[]		#testing variable used to make sure OSD_interpolator is producing reasonable snrs
  OSD_setup = OSD_data(allosds, kepids, periods_length, durations_length, grid)
  allosds = 0 # unload OSD table to save memory
  add_param_fixed(sim_param,"read_OSD_function",true)
  return OSD_setup
end

setup_OSD_interp(sim_param::SimParam; force_reread::Bool = false) = setup_OSD(sim_param, force_reread=force_reread)

function find_index_lower_bounding_point(grid::AbstractArray{T1,1}, x::T2; verbose::Bool = false) where {T1<:Real, T2<:Real}
    if verbose
        @assert issorted(grid)
        @assert length(grid)>=2
        println("# ", grid[1], " <= ", x, " <= ", grid[end])
        @assert grid[1] <= x <= grid[end]
    end
    idx = searchsortedlast(grid,x)
    if idx == 0
        idx = 1
    #elseif idx >= length(grid) # should never happen
    #    idx = length(grid)-1
    end
    return idx
end

function interp_OSD_from_table(kepid::Int64, period::T2, duration::T3; verbose::Bool = false) where {T2<:Real, T3<:Real}
  @assert eltype(OSD_setup.kepids) == Int64      # otherwise would need something like kepid = convert(eltype(OSD_setup.kepids),kepid)
  kepid_idx = searchsortedfirst(OSD_setup.kepids, kepid)
  if (kepid_idx > num_stars(OSD_setup)) || (OSD_setup.kepids[kepid_idx] != kepid) # if we don't find the kepid in allosds.jld, then we make a random one
     kepid_idx = rand(1:size(OSD_setup.allosds,1))
     if verbose
            println("# picked random kepid = ", OSD_setup.kepids[kepid_idx])
     end
  end
  idx_duration = find_index_lower_bounding_point(OSD_setup.grid[1], duration)
  idx_period   = find_index_lower_bounding_point(OSD_setup.grid[2], period)
  #z = view(OSD_setup.allosds,kepid_idx,idx_duration:(idx_duration+1),idx_period:(idx_period+1))     # use correct kepid index to extract 2D table from 3D OSD table

  #= value = z[1,1] * w_dur * w_per +
          z[2,1] * (1-w_dur) * w_per +
          z[1,2] * w_dur * (1-w_per) +
          z[2,2] * (1-w_dur) * (1-w_per)  =#

  if idx_duration < length(OSD_setup.grid[1]) && idx_period < length(OSD_setup.grid[2])
      w_dur = (duration-OSD_setup.grid[1][idx_duration]) / (OSD_setup.grid[1][idx_duration+1]-OSD_setup.grid[1][idx_duration])
      w_per = (period  -OSD_setup.grid[2][idx_period])   / (OSD_setup.grid[2][idx_period+1]  -OSD_setup.grid[2][idx_period])

      value = OSD_setup.allosds[kepid_idx,idx_duration,  idx_period  ] * w_dur * w_per +
          OSD_setup.allosds[kepid_idx,idx_duration+1,idx_period  ] * (1-w_dur) * w_per +
          OSD_setup.allosds[kepid_idx,idx_duration,  idx_period+1] * w_dur * (1-w_per) +
          OSD_setup.allosds[kepid_idx,idx_duration+1,idx_period+1] * (1-w_dur) * (1-w_per)
  elseif idx_period < length(OSD_setup.grid[2])
      w_per = (period  -OSD_setup.grid[2][idx_period])   / (OSD_setup.grid[2][idx_period+1]  -OSD_setup.grid[2][idx_period])

      value = (OSD_setup.allosds[kepid_idx,idx_duration,  idx_period+1] - OSD_setup.allosds[kepid_idx,idx_duration,  idx_period  ]) * w_per +
           OSD_setup.allosds[kepid_idx,idx_duration,  idx_period  ]
  elseif idx_duration < length(OSD_setup.grid[1])
      w_dur = (duration-OSD_setup.grid[1][idx_duration]) / (OSD_setup.grid[1][idx_duration+1]-OSD_setup.grid[1][idx_duration])

      value = (OSD_setup.allosds[kepid_idx,idx_duration+1,idx_period  ] - OSD_setup.allosds[kepid_idx,idx_duration,  idx_period  ]) * w_dur +
          OSD_setup.allosds[kepid_idx,idx_duration,idx_period  ]
  else
      value = OSD_setup.allosds[kepid_idx,idx_duration,  idx_period  ]
  end
end

# function interp_OSD_from_table(kepid::Int64, period::Real, duration::Real)
#   kepid = convert(Float64,kepid)
#   meskep = OSD_setup.kepids			#we need to find the index that this planet's kepid corresponds to in allosds.jld
#   kepid_index = findfirst(meskep, kepid)
#   if kepid_index == 0
#      kepid_index = rand(1:88807)		#if we don't find the kepid in allosds.jld, then we make a random one
#   end
#   olOSD = OSD_setup.allosds[kepid_index,:,:]    #use correct kepid index to extract 2D table from 3D OSD table
#   # olOSD = convert(Array{Float64,2},olOSD)
#   @time lint = Lininterp(olOSD, OSD_setup.grid)	#sets up linear interpolation object
#   osd = ApproXD.eval2D(lint, [duration*24,period])[1]	#interpolates osd
#   return osd
# end

# function cdpp_vs_osd(ratio::Float64, cuantos::Int64)
# #testing function that takes ratios of cdpp_snr/osd_snr and plots a histogram to make sure the results are reasonable.
#   global compareNoise
#   push!(compareNoise,ratio)
#   if length(compareNoise) == cuantos
#     PyPlot.plt[:hist](compareNoise,100)
#     println("MES median: ",median(compareNoise)," MES mean: ",mean(compareNoise), " Standard deviation: ",std(compareNoise))
#     cuantos = 100000000
#   end
#   return cuantos
# end

end  # module WindowFunction
