## ExoplanetsSysSim/examples/dr25_gaia_m/dr25_binrates_func.jl
## (c) 2019 Danley C. Hsu & Eric B. Ford
# Collection of functions specific to estimating DR25
#   planet candidate occurrence rates over a 2D period-radius grid

using ExoplanetsSysSim
using StatsFuns
using JLD
using CSV
using DataFrames
using Distributions

## simulation_parameters
macro isdefinedlocal(var) 
    quote 
        try 
            $(esc(var)) 
            true 
        catch err 
            isa(err, UndefVarError) ? false : rethrow(err) 
        end 
    end
end

function setup_sim_param_dr25binrates(args::Vector{String} = String[] )   # allow this to take a list of parameter (e.g., from command line)
    sim_param = ExoplanetsSysSim.SimParam()

    add_param_fixed(sim_param,"max_tranets_in_sys",7)
    add_param_fixed(sim_param,"generate_star",ExoplanetsSysSim.generate_star_dumb)
    add_param_fixed(sim_param,"generate_planetary_system", ExoplanetsSysSim.generate_planetary_system_uncorrelated_incl)
    add_param_fixed(sim_param,"generate_kepler_target",ExoplanetsSysSim.generate_kepler_target_from_table)
    add_param_fixed(sim_param,"star_table_setup",setup_star_table_dr25)
    add_param_fixed(sim_param,"stellar_catalog","q1q17_dr25_gaia_m.jld")
    add_param_fixed(sim_param,"osd_file","dr25m_osds.jld")
    add_param_fixed(sim_param,"generate_num_planets",generate_num_planets_binrates_uniform)
    add_param_fixed(sim_param,"generate_planet_mass_from_radius",ExoplanetsSysSim.generate_planet_mass_from_radius_powerlaw)
    add_param_fixed(sim_param,"vetting_efficiency",ExoplanetsSysSim.vetting_efficiency_none)  
    add_param_fixed(sim_param,"mr_power_index",2.0)
    add_param_fixed(sim_param,"mr_const",1.0)
    add_param_fixed(sim_param,"generate_period_and_sizes", generate_period_and_sizes_binrates_uniform)
    add_param_fixed(sim_param,"p_lim_full",[0.5, 1., 2., 4., 8., 16., 32., 64., 128., 256., 500.])
    add_param_fixed(sim_param,"r_lim_full",[0.25, 0.5, 0.75, 1., 1.25, 1.5, 1.75, 2., 2.5, 3., 4., 6., 8., 12., 16.]*ExoplanetsSysSim.earth_radius)
    #p_dim = length(p_lim_arr_num)-1
    #r_dim = length(r_lim_arr_num)-1
    #rate_tab_init = reshape(fill(1.0, p_dim*r_dim)*0.01,(r_dim,p_dim))
    #add_param_fixed(sim_param, "p_lim_arr", p_lim_arr_num)
    #add_param_fixed(sim_param, "r_lim_arr", r_lim_arr_num*ExoplanetsSysSim.earth_radius)
    #add_param_active(sim_param,"obs_par", rate_tab_init)
    add_param_fixed(sim_param,"generate_e_omega",ExoplanetsSysSim.generate_e_omega_rayleigh)
    add_param_fixed(sim_param,"sigma_hk",0.03)
    add_param_fixed(sim_param,"sigma_incl",2.0)   # degrees 
    add_param_fixed(sim_param,"calc_target_obs_sky_ave",ExoplanetsSysSim.calc_target_obs_sky_ave)
    add_param_fixed(sim_param,"calc_target_obs_single_obs",ExoplanetsSysSim.calc_target_obs_single_obs)
    add_param_fixed(sim_param,"transit_noise_model",ExoplanetsSysSim.transit_noise_model_diagonal)
    return sim_param
end

function set_test_param(sim_param_closure::SimParam)
    @eval(include(joinpath(pwd(),"param.in")))

    if @isdefinedlocal(stellar_catalog)
        @assert (typeof(stellar_catalog) == String)
        add_param_fixed(sim_param_closure,"stellar_catalog",stellar_catalog)
    end
    if @isdefinedlocal(koi_catalog)
        @assert (typeof(koi_catalog) == String)
        add_param_fixed(sim_param_closure,"koi_catalog",koi_catalog)
    end
    
    if @isdefinedlocal(num_targ_sim)
        @assert (typeof(num_targ_sim) == Int)
        add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",num_targ_sim)
    end

    if @isdefinedlocal(osd_file)
        @assert (typeof(osd_file) == String)
        add_param_fixed(sim_param_closure,"osd_file",osd_file)
    end
    
    @assert (typeof(p_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "p_lim_arr", p_bin_lim)

    @assert (typeof(r_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "r_lim_arr", r_bin_lim*ExoplanetsSysSim.earth_radius)

    p_dim = length(get_any(sim_param_closure, "p_lim_arr", Array{Float64,1}))-1
    r_dim = length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    n_bin = p_dim*r_dim

    if @isdefinedlocal(rate_init)
        if typeof(rate_init) <: Real
            @assert (rate_init >= 0.0)
            rate_init_list = fill(rate_init, n_bin)
        else
            rate_init_list = rate_init
        end
        
        @assert (ndims(rate_init_list) <= 2)
        if ndims(rate_init_list) == 1
            @assert (length(rate_init_list) == n_bin)
            rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        else
            @assert (size(rate_init_list) == (r_dim, p_dim))
            rate_tab_init = rate_init_list*0.01
        end
        add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    else
        rate_init_list = fill(1.0, n_bin)
        rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    end
    
    return sim_param_closure
end

function set_test_param_total(sim_param_closure::SimParam)
    @eval(include(joinpath(pwd(),"param.in")))

    if @isdefinedlocal(stellar_catalog)
        @assert (typeof(stellar_catalog) == String)
        add_param_fixed(sim_param_closure,"stellar_catalog",stellar_catalog)
    end
    if @isdefinedlocal(koi_catalog)
        @assert (typeof(koi_catalog) == String)
        add_param_fixed(sim_param_closure,"koi_catalog",koi_catalog)
    end
    
    if @isdefinedlocal(num_targ_sim)
        @assert (typeof(num_targ_sim) == Int)
        add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",num_targ_sim)
    end

    if @isdefinedlocal(osd_file)
        @assert (typeof(osd_file) == String)
        add_param_fixed(sim_param_closure,"osd_file",osd_file)
    end
    
    @assert (typeof(p_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "p_lim_arr", p_bin_lim)

    @assert (typeof(r_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "r_lim_arr", r_bin_lim*ExoplanetsSysSim.earth_radius)

    p_dim = length(get_any(sim_param_closure, "p_lim_arr", Array{Float64,1}))-1
    r_dim = length(get_any(sim_param_closure, "r_lim_arr", Array{Float64,1}))-1
    n_bin = p_dim*r_dim

    if @isdefinedlocal(rate_init)
        if typeof(rate_init) <: Real
            @assert (rate_init >= 0.0)
            rate_init_list = fill(rate_init, n_bin)
        else
            rate_init_list = rate_init
        end
        
        @assert (ndims(rate_init_list) <= 2)
        if ndims(rate_init_list) == 1
            @assert (length(rate_init_list) == n_bin)
            rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        else
            @assert (size(rate_init_list) == (r_dim, p_dim))
            rate_tab_init = rate_init_list*0.01
        end
        if r_dim > 1
            lamb_col = sum(rate_tab_init, dims=1)
            rate_tab_init = rate_tab_init ./ lamb_col
            rate_tab_init = vcat(lamb_col, rate_tab_init)
        end
        add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    else
        rate_init_list = fill(1.0, n_bin)
        rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        if r_dim > 1
            lamb_col = sum(rate_tab_init, dims=1)
            rate_tab_init = rate_tab_init ./ lamb_col
            rate_tab_init = vcat(lamb_col, rate_tab_init)
        end
        add_param_active(sim_param_closure, "obs_par", rate_tab_init)
    end
    
    if r_dim == 1
        add_param_fixed(sim_param_closure,"generate_period_and_sizes", generate_period_and_sizes_christiansen_single_rp)
    end
    
    return sim_param_closure
end

function set_test_param_ratio(sim_param_closure::SimParam)
    @eval(include(joinpath(pwd(),"param.in")))

    if @isdefinedlocal(stellar_catalog)
        @assert (typeof(stellar_catalog) == String)
        add_param_fixed(sim_param_closure,"stellar_catalog",stellar_catalog)
    end
    if @isdefinedlocal(koi_catalog)
        @assert (typeof(koi_catalog) == String)
        add_param_fixed(sim_param_closure,"koi_catalog",koi_catalog)
    end
    
    if @isdefinedlocal(num_targ_sim)
        @assert (typeof(num_targ_sim) == Int)
        add_param_fixed(sim_param_closure,"num_targets_sim_pass_one",num_targ_sim)
    end

    if @isdefinedlocal(osd_file)
        @assert (typeof(osd_file) == String)
        add_param_fixed(sim_param_closure,"osd_file",osd_file)
    end
    
    @assert (typeof(p_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "p_lim_arr", p_bin_lim)

    @assert (typeof(r_bin_lim) == Array{Float64,1})
    add_param_fixed(sim_param_closure, "r_lim_arr", r_bin_lim*ExoplanetsSysSim.earth_radius)

    p_dim = length(get_any(sim_param_closure, "p_lim_full", Array{Float64,1}))-1
    r_dim = length(get_any(sim_param_closure, "r_lim_full", Array{Float64,1}))-1
    n_bin = p_dim*r_dim

    if @isdefinedlocal(rate_init)
        if typeof(rate_init) <: Real
            @assert (rate_init >= 0.0)
            rate_init_list = fill(rate_init, n_bin)
        else
            rate_init_list = rate_init
        end
        
        @assert (ndims(rate_init_list) <= 2)
        if ndims(rate_init_list) == 1
            @assert (length(rate_init_list) == n_bin)
            rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        else
            @assert (size(rate_init_list) == (r_dim, p_dim))
            rate_tab_init = rate_init_list*0.01
        end
        add_param_fixed(sim_param_closure, "obs_par", rate_tab_init)
    else
        rate_init_list = fill(1.0, n_bin)
        rate_tab_init = reshape(rate_init_list*0.01, (r_dim, p_dim))
        add_param_fixed(sim_param_closure, "obs_par", rate_tab_init)
    end
    
    return sim_param_closure
end


## planetary_system
function draw_uniform_selfavoiding(n::Integer; lower_bound::Real=0.0, upper_bound=1.0, min_separation::Real = 0.05, return_sorted::Bool=false )
    @assert(n>=1)
    @assert(upper_bound>lower_bound)
    @assert(2*min_separation*n<upper_bound-lower_bound)
    list = rand(n)
    sorted_idx = collect(1:n)
    segment_length = upper_bound-lower_bound
    list[1] = lower_bound+segment_length*list[1]   # First draw is standard uniform
    segment_length -= min(upper_bound,list[1]+min_separation)-max(lower_bound,list[1]-min_separation) 
    for i in 2:n
        segment_length -= min(upper_bound,list[i-1]+min_separation)-max(lower_bound,list[i-1]-min_separation)   # Reduce length for future draws
        list[i] *= segment_length    # Draw over reduced range based on which segments need to be excluded
        list[i] += lower_bound
        j = 1
        while j<= i-1 # Checking for conflicts
            k = sorted_idx[j]     # Going from low to high
            if list[i]>list[k]-min_separation   # If too close, then bu
               list[i] += min(upper_bound,list[k]+min_separation)-max(lower_bound,list[k]-min_separation)
            else
                break
            end
            j += 1
        end
        for k in i:-1:(j+1)   # Keep larger values sorted
            sorted_idx[k]=sorted_idx[k-1]
        end
        sorted_idx[j] = i   # Save order for this draw
        #segment_length -= min(upper_bound,list[i]+min_separation)-max(lower_bound,list[i]-min_separation)   # Reduce length for future draws
   end 
   return return_sorted ? list[sorted_idx] : list
end

function generate_num_planets_binrates_uniform(s::Star, sim_param::SimParam)
  local max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys") # TODO SCI: Is 7 planets max per system OK, even when fitting across potentially 9 period bins?
  #local max_tranets_per_P::Int64 = 3  # Set maximum number of planets per period range as loose stability criteria and to prevent near-crossing orbits
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  local p_dim = length(limitP)-1
  local r_dim = length(get_any(sim_param, "r_lim_arr", Array{Float64,1}))-1
  sum_lambda = 0
  for i in 1:p_dim
      sum_lambda += ExoplanetsSysSim.generate_num_planets_poisson(sum(rate_tab[:,i]), convert(Int64, floor(3*log(limitP[i+1]/limitP[i])/log(2))))
  end
  #println("# lambda= ", sum_lambda) 
  return min(sum_lambda, max_tranets_in_sys)
end

function generate_num_planets_binrates_beta(s::Star, sim_param::SimParam)
  local max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys") # TODO SCI: Is 7 planets max per system OK, even when fitting across potentially 9 period bins?
  #local max_tranets_per_P::Int64 = 3  # Set maximum number of planets per period range as loose stability criteria and to prevent near-crossing orbits
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  local p_dim = length(limitP)-1
  local r_dim = length(get_any(sim_param, "r_lim_arr", Array{Float64,1}))-1
  sum_lambda = 0
  for i in 1:p_dim
      sum_lambda += ExoplanetsSysSim.generate_num_planets_poisson(rate_tab[1,i], convert(Int64, floor(3*log(limitP[i+1]/limitP[i])/log(2))))
  end
  #println("# lambda= ", sum_lambda) 
  return min(sum_lambda, max_tranets_in_sys)
end

function generate_num_planets_m_fgk_ratio(s::Star, sim_param::SimParam)
  local max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys") # TODO SCI: Is 7 planets max per system OK, even when fitting across potentially 9 period bins?
  #local max_tranets_per_P::Int64 = 3  # Set maximum number of planets per period range as loose stability criteria and to prevent near-crossing orbits
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  mfgk_ratio::Float64 = get_real(sim_param, "mfgk_ratio")  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_full", Array{Float64,1})
  local p_dim = length(limitP)-1
  local r_dim = length(get_any(sim_param, "r_lim_full", Array{Float64,1}))-1
  sum_lambda = 0
  for i in 1:p_dim
      sum_lambda += ExoplanetsSysSim.generate_num_planets_poisson(mfgk_ratio*sum(rate_tab[:,i]), convert(Int64, floor(3*log(limitP[i+1]/limitP[i])/log(2))))
  end
  #println("# lambda= ", sum_lambda) 
  return min(sum_lambda, max_tranets_in_sys)
end

function generate_period_and_sizes_binrates_uniform(s::Star, sim_param::SimParam; num_pl::Integer = 1)
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_arr", Array{Float64,1})
  local r_dim = length(limitRp)-1
  sepa_min = 0.05  # Minimum orbital separation in AU
  backup_sepa_factor_slightly_less_than_one = 0.95  
    
  @assert ((length(limitP)-1) == size(rate_tab, 2))
  @assert ((length(limitRp)-1) == size(rate_tab, 1))

  Plist = zeros(num_pl)
  Rplist = zeros(num_pl)
  rate_tab_1d = reshape(rate_tab,length(rate_tab))
  maxcuml = sum(rate_tab_1d)
  cuml = cumsum(rate_tab_1d/maxcuml)

  # We assume uniform sampling in log P and log Rp within each bin
  j_idx = ones(Int64, num_pl)
    
  for n in 1:num_pl
    rollp = Base.rand()
    idx = findfirst(x -> x > rollp, cuml)
    i_idx = (idx-1)%size(rate_tab,1)+1
    j_idx[n] = floor(Int64,(idx-1)//size(rate_tab,1))+1
    Rplist[n] = exp(Base.rand()*(log(limitRp[i_idx+1])-log(limitRp[i_idx]))+log(limitRp[i_idx]))
  end

  for j in 1:(length(limitP)-1)
      tmp_ind = findall(x -> x == j, j_idx)
      if length(tmp_ind) > 0
          redraw_att = 0
          invalid_config = true
          while invalid_config && redraw_att < 20
              n_range = length(tmp_ind)
              loga_min = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass))
              loga_min_ext = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass)+sepa_min)  # Used for determining minimum semimajor axis separation
              loga_max = log(ExoplanetsSysSim.semimajor_axis(limitP[j+1], s.mass))
              logsepa_min = min(loga_min_ext-loga_min, (loga_max-loga_min)/n_range/2*backup_sepa_factor_slightly_less_than_one)  # Prevents minimum separations too large
              tmp_logalist = draw_uniform_selfavoiding(n_range,min_separation=logsepa_min,lower_bound=loga_min,upper_bound=loga_max)
              tmp_Plist = exp.((3*tmp_logalist .- log(s.mass))/2)*ExoplanetsSysSim.day_in_year  # Convert from log a (in AU) back to P (in days)
              invalid_config = false
              redraw_att += 1
              for n in 1:n_range
                  if tmp_Plist[n] < limitP[j] || tmp_Plist[n] > limitP[j+1]
                      invalid_config = true
                  else
                      Plist[tmp_ind[n]] = tmp_Plist[n]
                  end
              end
          end
      end
  end
  return Plist, Rplist
end

function generate_period_and_sizes_binrates_dirichlet(s::Star, sim_param::SimParam; num_pl::Integer = 1)
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_arr", Array{Float64,1})
  local r_dim = length(limitRp)-1
  sepa_min = 0.05  # Minimum orbital separation in AU
  backup_sepa_factor_slightly_less_than_one = 0.95  
    
  @assert ((length(limitP)-1) == size(rate_tab, 2))
  @assert ((length(limitRp)-1) == (size(rate_tab, 1)-1))
    
  Plist = zeros(num_pl)
  Rplist = zeros(num_pl)
  maxcuml = sum(rate_tab[1,:])
  cuml = cumsum(rate_tab[1,:]/maxcuml)  

  # We assume uniform sampling in log P and log Rp within each bin
  j_idx = ones(Int64, num_pl)
    
  for n in 1:num_pl
    rollp = Base.rand()
    j_idx[n] = findfirst(x -> x > rollp, cuml)
  end
    
  for j in 1:(length(limitP)-1)
      tmp_ind = find(x -> x == j, j_idx)
      if length(tmp_ind) > 0
          redraw_att = 0
          invalid_config = true
          while invalid_config && redraw_att < 20
              n_range = length(tmp_ind)
              loga_min = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass))
              loga_min_ext = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass)+sepa_min)  # Used for determining minimum semimajor axis separation
              loga_max = log(ExoplanetsSysSim.semimajor_axis(limitP[j+1], s.mass))
              logsepa_min = min(loga_min_ext-loga_min, (loga_max-loga_min)/n_range/2*backup_sepa_factor_slightly_less_than_one)  # Prevents minimum separations too large
              tmp_logalist = draw_uniform_selfavoiding(n_range,min_separation=logsepa_min,lower_bound=loga_min,upper_bound=loga_max)
              tmp_Plist = exp.((3*tmp_logalist - log(s.mass))/2)*ExoplanetsSysSim.day_in_year  # Convert from log a (in AU) back to P (in days)
              rad_dist = Distributions.Categorical(rate_tab[((j-1)*(r_dim+1)+2):((j-1)*(r_dim+1)+(r_dim+1))]) # Distribution for fraction of times the next planet draw would be assigned to a given radius bin
              invalid_config = false
              redraw_att += 1
              for n in 1:n_range
                  if tmp_Plist[n] < limitP[j] || tmp_Plist[n] > limitP[j+1]
                      invalid_config = true
                  else
                      Plist[tmp_ind[n]] = tmp_Plist[n]
                  end
                  i_idx = rand(rad_dist)
                  Rplist[tmp_ind[n]] = exp(Base.rand()*(log(limitRp[i_idx+1])-log(limitRp[i_idx]))+log(limitRp[i_idx]))
              end
          end
      end
  end
  return Plist, Rplist
end

function generate_period_and_sizes_m_fgk_ratio(s::Star, sim_param::SimParam; num_pl::Integer = 1)
  rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
  #mfgk_ratio::Float64 = get_real(sim_param, "mfgk_ratio")
  
  limitP::Array{Float64,1} = get_any(sim_param, "p_lim_full", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(sim_param, "r_lim_full", Array{Float64,1})
  local r_dim = length(limitRp)-1
  sepa_min = 0.05  # Minimum orbital separation in AU
  backup_sepa_factor_slightly_less_than_one = 0.95  
    
  @assert ((length(limitP)-1) == size(rate_tab, 2))
  @assert ((length(limitRp)-1) == size(rate_tab, 1))

  Plist = zeros(num_pl)
  Rplist = zeros(num_pl)
  rate_tab_1d = reshape(rate_tab,length(rate_tab))
  maxcuml = sum(rate_tab_1d)
  cuml = cumsum(rate_tab_1d/maxcuml)

  # We assume uniform sampling in log P and log Rp within each bin
  j_idx = ones(Int64, num_pl)
    
  for n in 1:num_pl
    rollp = Base.rand()
    idx = findfirst(x -> x > rollp, cuml)
    i_idx = (idx-1)%size(rate_tab,1)+1
    j_idx[n] = floor(Int64,(idx-1)//size(rate_tab,1))+1
    Rplist[n] = exp(Base.rand()*(log(limitRp[i_idx+1])-log(limitRp[i_idx]))+log(limitRp[i_idx]))
  end

  for j in 1:(length(limitP)-1)
      tmp_ind = find(x -> x == j, j_idx)
      if length(tmp_ind) > 0
          redraw_att = 0
          invalid_config = true
          while invalid_config && redraw_att < 20
              n_range = length(tmp_ind)
              loga_min = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass))
              loga_min_ext = log(ExoplanetsSysSim.semimajor_axis(limitP[j], s.mass)+sepa_min)  # Used for determining minimum semimajor axis separation
              loga_max = log(ExoplanetsSysSim.semimajor_axis(limitP[j+1], s.mass))
              logsepa_min = min(loga_min_ext-loga_min, (loga_max-loga_min)/n_range/2*backup_sepa_factor_slightly_less_than_one)  # Prevents minimum separations too large
              tmp_logalist = draw_uniform_selfavoiding(n_range,min_separation=logsepa_min,lower_bound=loga_min,upper_bound=loga_max)
              tmp_Plist = exp.((3*tmp_logalist .- log(s.mass))/2)*ExoplanetsSysSim.day_in_year  # Convert from log a (in AU) back to P (in days)
              invalid_config = false
              redraw_att += 1
              for n in 1:n_range
                  if tmp_Plist[n] < limitP[j] || tmp_Plist[n] > limitP[j+1]
                      invalid_config = true
                  else
                      Plist[tmp_ind[n]] = tmp_Plist[n]
                  end
              end
          end
      end
  end
  return Plist, Rplist
end


## stellar_table
function setup_dr25(sim_param::SimParam; force_reread::Bool = false)
  #global df
  wf = WindowFunction.setup_window_function(sim_param)
  WindowFunction.setup_OSD_interp(sim_param) #read in osd files so they can be interpolated
  df = ExoplanetsSysSim.StellarTable.df
  if haskey(sim_param,"read_stellar_catalog") && !force_reread
     return df
     #return data
  end
  stellar_catalog_filename = convert(String,joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")), "data", convert(String,get(sim_param,"stellar_catalog","q1_q17_dr25_stellar.csv")) ) )
  df = setup_dr25(stellar_catalog_filename)
  add_param_fixed(sim_param,"read_stellar_catalog",true)
  add_param_fixed(sim_param,"num_kepler_targets",StellarTable.num_usable_in_star_table())
  if !haskey(sim_param.param,"num_targets_sim_pass_one")
      add_param_fixed(sim_param,"num_targets_sim_pass_one", StellarTable.num_usable_in_star_table())
  end
  StellarTable.set_star_table(df)
  return df  
end

function setup_dr25(filename::String; force_reread::Bool = false)
  #global df, usable
  df = ExoplanetsSysSim.StellarTable.df
  #usable = ExoplanetsSysSim.StellarTable.usable
  if occursin(r".jld2$",filename)
  try 
    data = load(filename)
    df = data["stellar_catalog"]
    #usable::Array{Int64,1} = data["stellar_catalog_usable"]
    Core.typeassert(df,DataFrame)
    StellarTable.set_star_table(df)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in jld2 format."))
  end

  else
  try 
    df = CSV.read(filename)
  catch
    error(string("# Failed to read stellar catalog >",filename,"< in ascii format."))
  end

  has_mass = .! (ismissing.(df[:mass]) .| ismissing.(df[:mass_err1]) .| ismissing.(df[:mass_err2]))
  has_radius = .! (ismissing.(df[:radius]) .| ismissing.(df[:radius_err1]) .| ismissing.(df[:radius_err2]))
  has_dens = .! (ismissing.(df[:dens]) .| ismissing.(df[:dens_err1]) .| ismissing.(df[:dens_err2]))
  has_cdpp = .! (ismissing.(df[:rrmscdpp01p5]) .| ismissing.(df[:rrmscdpp02p0]) .| ismissing.(df[:rrmscdpp02p5]) .| ismissing.(df[:rrmscdpp03p0]) .| ismissing.(df[:rrmscdpp03p5]) .| ismissing.(df[:rrmscdpp04p5]) .| ismissing.(df[:rrmscdpp05p0]) .| ismissing.(df[:rrmscdpp06p0]) .| ismissing.(df[:rrmscdpp07p5]) .| ismissing.(df[:rrmscdpp09p0]) .| ismissing.(df[:rrmscdpp10p5]) .| ismissing.(df[:rrmscdpp12p0]) .| ismissing.(df[:rrmscdpp12p5]) .| ismissing.(df[:rrmscdpp15p0]))
  has_ld = .! (ismissing.(df[:limbdark_coeff1]) .| ismissing.(df[:limbdark_coeff2]) .| ismissing.(df[:limbdark_coeff3]) .| ismissing.(df[:limbdark_coeff4]))
  has_rest = .! (ismissing.(df[:dataspan]) .| ismissing.(df[:dutycycle]))
  in_Q1Q12 = []
  obs_gt_5q = []
  for x in df[:st_quarters]
    subx = string(x)
    num_q_obs = length(matchall(r"1", subx))
    push!(obs_gt_5q, num_q_obs>5)
    subx = ("0"^(17-length(subx)))*subx
    indQ = search(subx, '1')
    if ((indQ < 1) | (indQ > 12)) || num_q_obs<=5 
      push!(in_Q1Q12, false)
    else
      push!(in_Q1Q12, true)
    end
  end
  is_FGK = []
  for x in 1:length(df[:teff])
    if ((df[x,:teff] > 4000.0) & (df[x,:teff] < 7000.0) & (df[x,:logg] > 4.0))
      push!(is_FGK, true)
    else
      push!(is_FGK, false)
    end
  end
  is_usable = has_radius .& is_FGK .& has_mass .& has_rest .& has_dens .& has_cdpp .& obs_gt_5q .& has_ld
  if contains(filename,"q1_q16_stellar.csv")
    is_usable = is_usable .& in_Q1Q12
  end
  # See options at: http://exoplanetarchive.ipac.caltech.edu/docs/API_keplerstellar_columns.html
  symbols_to_keep = [ :kepid, :mass, :mass_err1, :mass_err2, :radius, :radius_err1, :radius_err2, :dens, :dens_err1, :dens_err2, :rrmscdpp01p5, :rrmscdpp02p0, :rrmscdpp02p5, :rrmscdpp03p0, :rrmscdpp03p5, :rrmscdpp04p5, :rrmscdpp05p0, :rrmscdpp06p0, :rrmscdpp07p5, :rrmscdpp09p0, :rrmscdpp10p5, :rrmscdpp12p0, :rrmscdpp12p5, :rrmscdpp15p0, :cdppslplong, :cdppslpshrt, :dataspan, :dutycycle, :limbdark_coeff1, :limbdark_coeff2, :limbdark_coeff3, :limbdark_coeff4 ]
  delete!(df, [~(x in symbols_to_keep) for x in names(df)])    # delete columns that we won't be using anyway
  usable = findall(is_usable)
  df = df[usable, symbols_to_keep]
  tmp_df = DataFrame()    
  for col in names(df)
      tmp_df[col] = collect(skipmissing(df[col]))
  end
  df = tmp_df
  mast_df = CSV.read(convert(String,joinpath(abspath(joinpath(dirname(Base.find_package("ExoplanetsSysSim")),"..")), "data", "KeplerMAST_TargetProperties.csv")))
  delete!(mast_df, [~(x in [:kepid, :contam]) for x in names(mast_df)])
  df = join(df, mast_df, on=:kepid)
  StellarTable.set_star_table(df)
  end
    println("# Removing stars observed <5 quarters.")
    df[!,:wf_id] = map(x->ExoplanetsSysSim.WindowFunction.get_window_function_id(x,use_default_for_unknown=false),df[!,:kepid])
    obs_5q = df[!,:wf_id].!=-1
    df = df[obs_5q, names(df)]
    StellarTable.set_star_table(df)
  return df
end

setup_star_table_dr25(sim_param::SimParam; force_reread::Bool = false) = setup_dr25(sim_param, force_reread=force_reread)
setup_star_table_dr25(filename::String; force_reread::Bool = false) = setup_dr25(filename, force_reread=force_reread)


## summary_statistics
function calc_summary_stats_obs_binned_rates(cat_obs::KeplerObsCatalog, param::SimParam; trueobs_cat::Bool = false, obs_skyavg::Bool = false)
  ssd = Dict{String,Any}()
  cache = Dict{String,Any}()

  if !trueobs_cat
    ssd["num targets"] = get_int(param,"num_targets_sim_pass_one")
  else
    ssd["num targets"] = get_int(param,"num_kepler_targets")
  end

  max_tranets_in_sys = get_int(param,"max_tranets_in_sys")    # Demo that simulation parameters can specify how to evalute models, too
  @assert max_tranets_in_sys >= 1
  idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1}             # Find indices of systems with at least 1 tranet = potentially detectable transiting planet

  # Count total number of tranets and compile indices for N-tranet systems
  num_tranets = 0
  idx_n_tranets = Vector{Int64}[ Int64[] for m = 1:max_tranets_in_sys]
  for n in 1:max_tranets_in_sys-1
    idx_n_tranets[n] = findall(x::KeplerTargetObs-> length(x.obs) == n, cat_obs.target[idx_tranets] )
    num_tranets += n*length(idx_n_tranets[n])
  end
  idx_n_tranets[max_tranets_in_sys] = findall(x::KeplerTargetObs-> length(x.obs) >= max_tranets_in_sys, cat_obs.target[idx_tranets] )

  num_tranets += max_tranets_in_sys*length(idx_n_tranets[max_tranets_in_sys])  # WARNING: this means we need to ignore planets w/ indices > max_tranets_in_sys
  if ( length( findall(x::KeplerTargetObs-> length(x.obs) > max_tranets_in_sys, cat_obs.target[idx_tranets] ) ) > 0)   # Make sure max_tranets_in_sys is at least big enough for observed systems
    warn("Observational data has more transiting planets in one systems than max_tranets_in_sys allows.")
  end
  num_tranets  = convert(Int64,num_tranets)            # TODO OPT: Figure out why isn't this already an Int.  I may be doing something that prevents some optimizations

  num_sys_tranets = zeros(max_tranets_in_sys)                           # Since observed data, don't need to calculate probabilities.
  for n in 1:max_tranets_in_sys                                         # Make histogram of N-tranet systems
    num_sys_tranets[n] = length(idx_n_tranets[n])
  end
  ssd["num_sys_tranets"] = num_sys_tranets
  ssd["planets detected"] = num_tranets 

  period_list = zeros(num_tranets)
  weight_list = zeros(num_tranets)
  radius_list = zeros(num_tranets)

  n = 1    # tranet id

  if !trueobs_cat
    for i in idx_tranets
      ld = ExoplanetsSysSim.LimbDarkeningParam4thOrder(ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff1), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff2), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff3), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff4) )
      flux_ratio = (1.0+ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id, :contam))/1.0 # WARNING: Assumes flux = 1
      #Rstar = trueobs_cat ? cat_obs.target[i].star.radius : ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id, :radius) 
      for j in 1:num_planets(cat_obs.target[i])
        period_list[n] = cat_obs.target[i].obs[j].period
        if obs_skyavg  
          weight_list[n] = min(ExoplanetsSysSim.prob_detect(cat_obs.target[i].prob_detect,j), 1.0) # CHECK WHAT THIS DOES
        else
          weight_list[n] = 1.0
        end
        radius_ratio = ExoplanetsSysSim.ratio_from_depth(cat_obs.target[i].obs[j].depth*flux_ratio, ld)
        radius_list[n] = radius_ratio*ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id, :radius)
        #radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
        #radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*Rstar
        n = n+1
      end
    end
  else
    for i in idx_tranets
      ld = ExoplanetsSysSim.LimbDarkeningParam4thOrder(ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff1), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff2), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff3), ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id,:limbdark_coeff4) )
      flux_ratio = (1.0+ExoplanetsSysSim.StellarTable.star_table(cat_obs.target[i].star.id, :contam))/1.0 # WARNING: Assumes flux = 1  
      for j in 1:num_planets(cat_obs.target[i])
        period_list[n] = cat_obs.target[i].obs[j].period
        weight_list[n] = 1.0
        radius_ratio = ExoplanetsSysSim.ratio_from_depth(cat_obs.target[i].obs[j].depth*flux_ratio, ld)
        radius_list[n] = radius_ratio*cat_obs.target[i].star.radius            
        #radius_list[n] = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
        n = n+1
      end
    end
  end

  #ssd["period_list"] = period_list
  ssd["weight_list"] = weight_list
  #ssd["radius_list"] = radius_list
  
  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
  np_bin_idx = 1
  bin_match_list = fill(fill(0,0),(length(limitP)-1)*(length(limitRp)-1))
  
  for i in 1:(length(limitP)-1)
    P_match = findall(x -> ((x > limitP[i]) && (x < limitP[i+1])), period_list)
    for j in 1:(length(limitRp)-1)
      R_match = findall(x -> ((x > limitRp[j]) && (x < limitRp[j+1])), radius_list)
      
      bin_match = intersect(P_match, R_match)  
      bin_match_list[np_bin_idx] = bin_match
      np_bin[np_bin_idx] = sum(weight_list[bin_match])
      np_bin_idx += 1
    end
  end

  cache["bin_match_list"] = bin_match_list
  #ssd["planets detected"] = sum(np_bin)
  ssd["planets table"] = np_bin

  return CatalogSummaryStatistics(ssd, cache)
end


## abc_distance
function calc_distance_vector_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, pass::Int64, sim_param::SimParam ; verbose::Bool = false)
  p_dim = length(get_any(sim_param, "p_lim_arr", Array{Float64,1}))-1
  r_dim = length(get_any(sim_param, "r_lim_arr", Array{Float64,1}))-1
  #rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
    
  d = Array{Float64}(undef,0)
  if pass == 1
    if verbose
      println("# Summary 1, pass 1: ",summary1)
      println("# Summary 2, pass 1: ",summary2)
    end
    d = zeros(3)

    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]
    np_bin = zeros(length(np1))
    num_detect_sim = zeros(length(np1))

    ### Bernoulli distance
    bin_match_list = summary2.cache["bin_match_list"]
    @assert length(bin_match_list) == length(np1) 
    np2 = zeros(Int64,length(np1))
    ###  

      for n in 1:length(np1)
        #np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
        #np_bin[n] = dist_L2_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
        #np_bin[n] = distance_poisson_draw(np2[n]/summary2.stat["num targets"]*summary1.stat["num targets"], convert(Int64, np1[n]))
        np_bin[n], num_detect_sim[n] = distance_sum_of_bernoulli_draws(floor(Int64,np1[n]),summary1.stat["num targets"], summary2.stat["weight_list"], summary2.stat["num targets"], bin_match_list[n])

      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end

      #d[1] = maximum(np_bin)
      #d[1] = sum(np_bin)
      np1_ratio = np1 ./ summary1.stat["num targets"]
      np2_ratio = num_detect_sim ./ summary1.stat["num targets"]
      d[1] = distance_canberra(np1_ratio, np2_ratio)# + distance_cosine(np1_ratio, np2_ratio)
      
      #println("Total rate: ", rate_tab[1,1], " / Distance (radii): ", d[1], " / Sim. cat. ratio = ", sum(num_detect_sim[1:r_dim])/summary2.stat["num_targets"], " / Obs. cat. ratio = ", sum(np1[1:r_dim])/summary1.stat["num targets"], " / Distance (total): ", dist_L2_abs(sum(num_detect_sim[1:r_dim])/summary2.stat["num targets"], sum(np1[1:r_dim])/summary1.stat["num targets"])*r_dim)
      
       # for j in 1:p_dim
       #     d[1] += dist_L2_abs(sum(num_detect_sim[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])/summary1.stat["num targets"], sum(np1[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])/summary1.stat["num targets"])*r_dim
       # end
    else
    println("# calc_distance_vector_demo doesn't know what to do for pass= ", pass)
  end
  return d
end

#=
function calc_distance_vector_relative_binned(summary1::CatalogSummaryStatistics, summary2::CatalogSummaryStatistics, pass::Int64, sim_param::SimParam ; verbose::Bool = false)
  p_dim = length(get_any(sim_param, "p_lim_arr", Array{Float64,1}))-1
  r_dim = length(get_any(sim_param, "r_lim_arr", Array{Float64,1}))-1
  #rate_tab::Array{Float64,2} = get_any(sim_param, "obs_par", Array{Float64,2})
    
  d = Array{Float64}(0)
  if pass == 1
    if verbose
      println("# Summary 1, pass 1: ",summary1)
      println("# Summary 2, pass 1: ",summary2)
    end
    d = zeros(3)

    np1 = haskey(summary1.stat,"planets table") ? summary1.stat["planets table"] : summary1.stat["expected planets table"]
    np2 = haskey(summary2.stat,"planets table") ? summary2.stat["planets table"] : summary2.stat["expected planets table"]
    np_bin = zeros(length(np1))
    num_detect_sim = zeros(length(np1))

    ### Bernoulli distance
    bin_match_list = summary2.cache["bin_match_list"]
    @assert length(bin_match_list) == length(np1) 
    np2 = zeros(Int64,length(np1))
    ###  

      for n in 1:length(np1)
        #np_bin[n] = dist_L1_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
        #np_bin[n] = dist_L2_abs(np1[n]/summary1.stat["num targets"], np2[n]/summary2.stat["num targets"])
        #np_bin[n] = distance_poisson_draw(np2[n]/summary2.stat["num targets"]*summary1.stat["num targets"], convert(Int64, np1[n]))
        np_bin[n], num_detect_sim[n] = distance_sum_of_bernoulli_draws(floor(Int64,np1[n]),summary1.stat["num targets"], summary2.stat["weight_list"], summary2.stat["num targets"], bin_match_list[n])

      #println("True # [Bin ", n,"] = ",np1[n],", Expected # [Bin ", n,"] = ",np2[n])
    end
      #d[1] = maximum(np_bin)
      d[1] = 0.0
      np1_col = 0
      np2_col = 0
      #println("Total rate: ", rate_tab[1,1], " / Distance (radii): ", d[1], " / Sim. cat. ratio = ", sum(num_detect_sim[1:r_dim])/summary2.stat["num_targets"], " / Obs. cat. ratio = ", sum(np1[1:r_dim])/summary1.stat["num targets"], " / Distance (total): ", dist_L2_abs(sum(num_detect_sim[1:r_dim])/summary2.stat["num targets"], sum(np1[1:r_dim])/summary1.stat["num targets"])*r_dim)
      
      for j in 1:p_dim
          d[1] += dist_L2_abs(sum(num_detect_sim[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])/summary1.stat["num targets"], sum(np1[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])/summary1.stat["num targets"])*r_dim
          np1_col = sum(np1[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])
          np2_col = sum(num_detect_sim[(j-1)*r_dim+1:(j-1)*r_dim+r_dim])
          for i in 1:r_dim
              d[1] += dist_L2_abs(num_detect_sim[(j-1)*r_dim+i]/np2_col, np1[(j-1)*r_dim+i]/np1_col)
          end
      end
    else
    println("# calc_distance_vector_demo doesn't know what to do for pass= ", pass)
  end
  return d
end
=#

## eval_model
# function test_dr25binrates()
#   global sim_param_closure = setup_sim_param_dr25binrates()
#   cat_phys = generate_kepler_physical_catalog(sim_param_closure)
#   cat_obs = observe_kepler_targets_single_obs(cat_phys,sim_param_closure)
#   global summary_stat_ref_closure = calc_summary_stats_obs_demo(cat_obs,sim_param_closure)
#   global cat_phys_try_closure  = generate_christiansen_catalog(sim_param_closure)
#   global cat_obs_try_closure  = observe_kepler_targets_sky_avg(cat_phys_try_closure,sim_param_closure)
#   global summary_stat_try_closure  = calc_summary_stats_sim_pass_one_demo(cat_obs_try_closure,cat_phys_try_closure,sim_param_closure)
#   summary_stat_try_closure   = calc_summary_stats_sim_pass_two_demo(cat_obs_try_closure,cat_phys_try_closure,summary_stat_try_closure,sim_param_closure)
#   param_guess = make_vector_of_sim_param(sim_xparam_closure)
#   evaluate_model_scalar_ret( param_guess)
# end


## inverse_detection & simple bayesian
function inv_det(cat_obs::KeplerObsCatalog, param::SimParam)
    num_targ = ExoplanetsSysSim.StellarTable.num_usable_in_star_table()

    limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
    limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

    println("------------------------------")
    cnt_bin, np_bin = cnt_np_bin(cat_obs, param)
    println("------------------------------")

    println("Inverse Detection Rates:")
    for i in 1:(length(limitP)-1)
        for j in 1:(length(limitRp)-1)
            rate_f = np_bin[(i-1)*(length(limitRp)-1) + j]/num_targ*100.
            if cnt_bin[(i-1)*(length(limitRp)-1) + j] > 0.
                println(rate_f, 
                        " +/- ", rate_f/sqrt(cnt_bin[(i-1)*(length(limitRp)-1) + j]), " %")
            else
                println(rate_f, 
                        " +/- N/A %")
            end
        end
    end
    println()
end

function simp_bayes(cat_obs::KeplerObsCatalog, param::SimParam)
    num_targ = ExoplanetsSysSim.StellarTable.num_usable_in_star_table()

    limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
    limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

    println("------------------------------")
    cnt_bin, np_bin = cnt_np_bin(cat_obs, param)
    println("------------------------------")
    ess_bin = stellar_ess(param)
    println("------------------------------")

    println("Simple Bayesian Rates:")
    for i in 1:(length(limitP)-1)
        for j in 1:(length(limitRp)-1)
            rate_f = (1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j])/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])*100.
            up_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.8413)*100.
            low_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.1587)*100.
            println(rate_f, 
                    " + ", up_quant - rate_f,
                    " - ", rate_f - low_quant, " %")
        end
    end
    println()
end

function inv_det_simp_bayes(cat_obs::KeplerObsCatalog, param::SimParam)
    num_targ = ExoplanetsSysSim.StellarTable.num_usable_in_star_table()

    limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
    limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

    println("------------------------------")
    cnt_bin, np_bin = cnt_np_bin(cat_obs, param)
    println("------------------------------")
    ess_bin = stellar_ess(param)
    println("------------------------------")

    println("Inverse Detection Rates:")
    for i in 1:(length(limitP)-1)
        for j in 1:(length(limitRp)-1)
            rate_f = np_bin[(i-1)*(length(limitRp)-1) + j]/num_targ*100.
            if cnt_bin[(i-1)*(length(limitRp)-1) + j] > 0.
                println(rate_f, 
                        " +/- ", rate_f/sqrt(cnt_bin[(i-1)*(length(limitRp)-1) + j]), " %")
            else
                println(rate_f, 
                        " +/- N/A %")
            end
        end
    end

    println()
    println("Simple Bayesian Rates:")
    for i in 1:(length(limitP)-1)
        for j in 1:(length(limitRp)-1)
            rate_f = (1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j])/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])*100.
            up_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.8413)*100.
            low_quant = quantile(Gamma(1.0+cnt_bin[(i-1)*(length(limitRp)-1) + j], 1.0/(1.0+ess_bin[(i-1)*(length(limitRp)-1) + j])), 0.1587)*100.
            println(rate_f, 
                    " + ", up_quant - rate_f,
                    " - ", rate_f - low_quant, " %")
        end
    end
    println()
end

## cnt_bin & np_bin (inverse detection & simple bayesian)
function cnt_np_bin(cat_obs::KeplerObsCatalog, param::SimParam, verbose::Bool = true)
    num_targ = ExoplanetsSysSim.StellarTable.num_usable_in_star_table()
    idx_tranets = findall(x::KeplerTargetObs-> length(x.obs) > 0, cat_obs.target)::Array{Int64,1} 

    limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
    limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

    np_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
    cnt_bin = zeros((length(limitP)-1) * (length(limitRp)-1))
    pl_idx = 1

    println("Calculating completeness for each planet...")
    for i in idx_tranets
        for j in 1:num_planets(cat_obs.target[i])
            pper = cat_obs.target[i].obs[j].period
            prad = sqrt(cat_obs.target[i].obs[j].depth)*cat_obs.target[i].star.radius
            
            pbin = findfirst(x -> ((pper > limitP[x]) && (pper < limitP[x+1])), collect(1:(length(limitP)-1)))
            rbin = findfirst(x -> ((prad > limitRp[x]) && (prad < limitRp[x+1])), collect(1:(length(limitRp)-1)))
            if (pbin > 0 && rbin > 0)
                cnt_bin[(pbin-1)*(length(limitRp)-1) + rbin] += 1
                pgeo = ExoplanetsSysSim.calc_transit_prob_single_planet_approx(pper, cat_obs.target[i].star.radius, cat_obs.target[i].star.mass)
	        pdet = 0.0
	        for star_id in 1:num_targ
                    ld = ExoplanetsSysSim.LimbDarkeningParam4thOrder(ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff1), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff2), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff3), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff4) )
	            star = SingleStar(ExoplanetsSysSim.StellarTable.star_table(star_id,:radius),ExoplanetsSysSim.StellarTable.star_table(star_id,:mass),1.0, ld, star_id)
                    cdpp_arr = ExoplanetsSysSim.make_cdpp_array_empty(star_id)#(1.0e-6*sqrt(1.0 / 24.0 / ExoplanetsSysSim.LC_duration)) .* [ExoplanetsSysSim.StellarTable.star_table(star_id, :rrmscdpp01p5)*sqrt(1.5), ExoplanetsSysSim.StellarTable.star_table(star_id, :rrmscdpp02p0)*sqrt(2.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp02p5)*sqrt(2.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp03p0)*sqrt(3.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp03p5)*sqrt(3.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp04p5)*sqrt(4.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp05p0)*sqrt(5.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp06p0)*sqrt(6.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp07p5)*sqrt(7.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp09p0)*sqrt(9.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp10p5)*sqrt(10.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp12p0)*sqrt(12.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp12p5)*sqrt(12.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp15p0)*sqrt(15.)]
                    #cdpp = 1.0e-6 * ExoplanetsSysSim.StellarTable.star_table(star_id, :rrmscdpp04p5) * sqrt(4.5/24.0 / ExoplanetsSysSim.LC_duration )
	            contam = 0.0
	            data_span = ExoplanetsSysSim.StellarTable.star_table(star_id, :dataspan)
	            duty_cycle = ExoplanetsSysSim.StellarTable.star_table(star_id, :dutycycle)
	            pl_arr = Array{Planet}(undef,1)
	            orbit_arr = Array{Orbit}(undef,1)
                    incl = acos(Base.rand()*star.radius*ExoplanetsSysSim.rsol_in_au/ExoplanetsSysSim.semimajor_axis(pper, star.mass))
	            orbit_arr[1] = Orbit(pper, 0., incl, 0., 0., Base.rand()*2.0*pi)
	            pl_arr[1] = Planet(prad, 1.0e-6)
                    if ExoplanetsSysSim.StellarTable.star_table_has_key(:wf_id)
                        wf_id = ExoplanetsSysSim.StellarTable.star_table(star_id,:wf_id)
                    else
                        wf_id = ExoplanetsSysSim.WindowFunction.get_window_function_id(ExoplanetsSysSim.StellarTable.star_table(star_id,:kepid))
                    end
	            kep_targ = KeplerTarget([PlanetarySystem(star, pl_arr, orbit_arr)], repeat(cdpp_arr, outer=[1,1]),contam,data_span,duty_cycle,wf_id)
                    
	            duration = ExoplanetsSysSim.calc_transit_duration(kep_targ,1,1) 
	            if duration <= 0.
	                continue
	            end
	            ntr = ExoplanetsSysSim.calc_expected_num_transits(kep_targ, 1, 1, param)
	            depth = ExoplanetsSysSim.calc_transit_depth(kep_targ,1,1)
                    cdpp = ExoplanetsSysSim.interpolate_cdpp_to_duration_lookup_cdpp(kep_targ, duration)
                    #cdpp = ExoplanetsSysSim.interpolate_cdpp_to_duration(kep_targ, duration)
                    snr = ExoplanetsSysSim.calc_snr_if_transit_cdpp(kep_targ, depth, duration, cdpp, param, num_transit=ntr)
	            pdet += ExoplanetsSysSim.calc_prob_detect_if_transit(kep_targ, snr, pper, duration, param, num_transit=ntr)
	        end
                np_bin[(pbin-1)*(length(limitRp)-1) + rbin] += 1.0/pgeo/(pdet/num_targ)
                if verbose
	            println("Planet ",pl_idx," => Bin ", (pbin-1)*(length(limitRp)-1) + rbin, ", C = ", 1.0/pgeo/(pdet/num_targ))
                end
	        pl_idx += 1
            end
        end
    end
    return cnt_bin, np_bin
end

## stellar catalog ess (simple bayesian)
function stellar_ess(param::SimParam, verbose::Bool = true)
  num_realiz = 100
  num_targ = ExoplanetsSysSim.StellarTable.num_usable_in_star_table()

  limitP::Array{Float64,1} = get_any(param, "p_lim_arr", Array{Float64,1})
  limitRp::Array{Float64,1} = get_any(param, "r_lim_arr", Array{Float64,1})

  ess_bin = zeros((length(limitP)-1) * (length(limitRp)-1))

  println(string("Stellar ESS calculation beginning..."))
  for star_id in 1:num_targ
    ld = ExoplanetsSysSim.LimbDarkeningParam4thOrder(ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff1), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff2), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff3), ExoplanetsSysSim.StellarTable.star_table(star_id,:limbdark_coeff4) )    
    star = SingleStar(ExoplanetsSysSim.StellarTable.star_table(star_id,:radius),ExoplanetsSysSim.StellarTable.star_table(star_id,:mass),1.0, ld, star_id)
    cdpp_arr = ExoplanetsSysSim.make_cdpp_array_empty(star_id)#(1.0e-6*sqrt(1.0/24.0/ExoplanetsSysSim.LC_duration)) .* [ExoplanetsSysSim.StellarTable.star_table(star_id, :rrmscdpp01p5)*sqrt(1.5), ExoplanetsSysSim.StellarTable.star_table(star_id, :rrmscdpp02p0)*sqrt(2.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp02p5)*sqrt(2.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp03p0)*sqrt(3.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp03p5)*sqrt(3.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp04p5)*sqrt(4.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp05p0)*sqrt(5.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp06p0)*sqrt(6.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp07p5)*sqrt(7.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp09p0)*sqrt(9.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp10p5)*sqrt(10.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp12p0)*sqrt(12.), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp12p5)*sqrt(12.5), ExoplanetsSysSim.StellarTable.star_table(star_id,:rrmscdpp15p0)*sqrt(15.)]
    contam = 0.0
    data_span = ExoplanetsSysSim.StellarTable.star_table(star_id, :dataspan)
    duty_cycle = ExoplanetsSysSim.StellarTable.star_table(star_id, :dutycycle)
    if ExoplanetsSysSim.StellarTable.star_table_has_key(:wf_id)
        wf_id = ExoplanetsSysSim.StellarTable.star_table(star_id,:wf_id)
    else
        wf_id = ExoplanetsSysSim.WindowFunction.get_window_function_id(ExoplanetsSysSim.StellarTable.star_table(star_id,:kepid))
    end
      
    for i_idx in 1:(length(limitP)-1)
      for j_idx in 1:(length(limitRp)-1)
        temp_bin = 0.0
        for n_test in 1:num_realiz
          pper = exp(Base.rand()*(log(limitP[i_idx+1])-log(limitP[i_idx]))+log(limitP[i_idx]))
	  prad = exp(Base.rand()*(log(limitRp[j_idx+1])-log(limitRp[j_idx]))+log(limitRp[j_idx]))
    
          pgeo = ExoplanetsSysSim.calc_transit_prob_single_planet_approx(pper, star.radius, star.mass)
	  pdet = 0.0
	
	  pl_arr = Array{Planet}(undef,1)
	  orbit_arr = Array{Orbit}(undef,1)
          incl = acos(Base.rand()*star.radius*ExoplanetsSysSim.rsol_in_au/ExoplanetsSysSim.semimajor_axis(pper, star.mass))
	  orbit_arr[1] = Orbit(pper, 0., incl, 0., 0., Base.rand()*2.0*pi)
	  pl_arr[1] = Planet(prad, 1.0e-6)  
	  kep_targ = KeplerTarget([PlanetarySystem(star, pl_arr, orbit_arr)], repeat(cdpp_arr, outer=[1,1]),contam,data_span,duty_cycle,wf_id)

	  duration = ExoplanetsSysSim.calc_transit_duration(kep_targ,1,1) 
	  if duration <= 0.
	    continue
	  end
	  ntr = ExoplanetsSysSim.calc_expected_num_transits(kep_targ, 1, 1, param)
	  depth = ExoplanetsSysSim.calc_transit_depth(kep_targ,1,1)
          # Apply correction to snr if grazing transit
          size_ratio = kep_targ.sys[1].planet[1].radius/kep_targ.sys[1].star.radius
          b = ExoplanetsSysSim.calc_impact_parameter(kep_targ.sys[1],1)
          snr_correction = ExoplanetsSysSim.calc_depth_correction_for_grazing_transit(b,size_ratio)
          depth *= snr_correction

          #cdpp = ExoplanetsSysSim.interpolate_cdpp_to_duration(kep_targ, duration)
          cdpp = ExoplanetsSysSim.interpolate_cdpp_to_duration_lookup_cdpp(kep_targ, duration)
          snr = ExoplanetsSysSim.calc_snr_if_transit_cdpp(kep_targ, depth, duration, cdpp, param, num_transit=ntr)
          #kepid = ExoplanetsSysSim.StellarTable.star_table(kep_targ.sys[1].star.id, :kepid)
          #osd_duration = ExoplanetsSysSim.get_legal_durations(pper,duration)	#tests if durations are included in Kepler's observations for a certain planet period. If not, returns nearest possible duration
          #osd = ExoplanetsSysSim.WindowFunction.interp_OSD_from_table(kepid, pper, osd_duration)
          #if osd_duration > duration				#use a correcting factor if this duration is lower than the minimum searched for this period. 
	  #   osd = osd*osd_duration/duration
          #end
          #snr = ExoplanetsSysSim.calc_snr_if_transit(kep_targ, depth, duration, osd, sim_param, num_transit=ntr)
	  pdet = ExoplanetsSysSim.calc_prob_detect_if_transit(kep_targ, snr, pper, duration, param, num_transit=ntr)
	
	  temp_bin += (pgeo*pdet)          
        end
        ess_bin[(i_idx-1)*(length(limitRp)-1) + j_idx] += temp_bin/num_realiz
      end
    end
    if verbose && rem(star_id, 10^convert(Int,floor(log10(num_targ)))) == 0.
      println(string("Star #", star_id, " finished"))
    end
  end

  if verbose
      println("")
      for i in 1:(length(limitP)-1)
          for j in 1:(length(limitRp)-1)
              println("Period limits: ", limitP[i:i+1], " / Radius limits: ", limitRp[j:j+1]/ExoplanetsSysSim.earth_radius, " / Stellar ESS = ", ess_bin[(i-1)*(length(limitRp)-1) + j])
          end
      end
  end
  return ess_bin
end
