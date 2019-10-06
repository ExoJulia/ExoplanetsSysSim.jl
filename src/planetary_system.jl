#using Distributions
#include("constants.jl")
#include("orbit.jl")
#include("planet.jl")

#if !@isdefined PlanetarySystemAbstract
  @compat abstract type PlanetarySystemAbstract end

  struct PlanetarySystem{StarT<:StarAbstract} <: PlanetarySystemAbstract
    star::StarT
    planet::Vector{Planet}
    orbit::Vector{Orbit}

      # TODO DETAIL: Setup inner constructor to enforce equal number of planets & orbits
      function PlanetarySystem{StarT}(s::StarT, p::AbstractVector{Planet}, o::AbstractVector{Orbit}) where {StarT<:StarAbstract}
        @assert(length(p)==length(o)) # else error(string("Number of planets must match number of orbits: Np= ",length(p)," No= ",length(o)))
        new(s,p,o)
      end
  end

  const PlanetarySystemSingleStar = PlanetarySystem{SingleStar}

#end

function PlanetarySystem(s::StarT) where {StarT<:StarAbstract}
     PlanetarySystem(s,Vector{Planet}(undef,0),Vector{Orbit}(undef,0))  # Constructor for a Planetary System with no planets
end

function PlanetarySystem(s::StarT, p::Planet, o::Orbit) where {StarT<:StarAbstract}
   PlanetarySystem(s,[p],[o])  # Constructor for a single Planet System
end

function PlanetarySystem(s::StarT, p::AbstractVector{Planet}, o::AbstractVector{Orbit}) where {StarT<:StarAbstract}
   PlanetarySystem{StarT}(s,p,o)  # Constructor for a single Planet System
end

function PlanetarySystem(ps::PlanetarySystem{StarT}, keep::AbstractVector{Int64})  where {StarT<:StarAbstract} # Why doesn't this work?
   PlanetarySystem{StarT}(ps.star,ps.planet[keep],ps.orbit[keep])
end

function star( ps::PlanetarySystem{StarT} )::StarT where {StarT<:StarAbstract}
  return ps.star
end

function planets( ps::PlanetarySystem{StarT} )::Vector{Planet} where {StarT<:StarAbstract}
  return ps.planet
end

function orbits( ps::PlanetarySystem{StarT} )::Vector{Orbit} where {StarT<:StarAbstract}
  return ps.orbit
end

#function PlanetarySystemSingleStar(ps::PlanetarySystemSingleStar, keep::Vector{Int64})
#function PlanetarySystem{StarT<:StarAbstract}(ps::PlanetarySystem{StarT}, keep::Vector{Int64})
#   PlanetarySystem(ps.star,ps.planet[keep],ps.orbit[keep])
#end

flux(ps::PlanetarySystem{StarT}) where {StarT<:StarAbstract} = flux(star(ps))
#flux(ps::PlanetarySystem{Star}) = flux(ps.star)
#flux(ps::PlanetarySystem{BinaryStar}) = flux(ps.star)
#flux(ps::PlanetarySystem{MultipleStar}) = flux(ps.star)

function num_planets(s::PlanetarySystem{StarT}) where {StarT<:StarAbstract}
  @assert( length(planets(s)) == length(orbits(s)) )    # TODO OPT: Deactivate inner assert's like this for speed once tested
  return length(planets(s))
end


function calc_hill_sphere(a::Float64, mu::Float64)
    return a*(mu/3)^(1//3)
end

function calc_mutual_hill_radii(ps::PlanetarySystem{StarT}, pl1::Int64, pl2::Int64) where StarT <: StarAbstract
    mu = (ps.planet[pl1].mass + ps.planet[pl2].mass)/ps.star.mass
    a = 0.5*(ps.orbit[pl1].a + ps.orbit[pl2].a)
    return calc_hill_sphere(a, mu)
end

function test_stability_circular(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, star_mass::Float64, sim_param::SimParam)
    @assert length(P) == length(mass)
    min_num_mutual_hill_radii = get_real(sim_param, "num_mutual_hill_radii")
    found_instability = false
    order = sortperm(P)
    a2 = semimajor_axis(P[order[1]], star_mass)
    for pl in 1:(length(P)-1)
        a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
        a2 = semimajor_axis(P[order[pl+1]], star_mass)
        a = 0.5*(a1+a2)
        mu = (mass[order[pl]] + mass[order[pl+1]])/star_mass
        mutual_hill_radius = calc_hill_sphere(a, mu)
        if a2-a1  < min_num_mutual_hill_radii*mutual_hill_radius
            found_instability = true
            break
        end
    end # loop over neighboring planet pairs within cluster
    return !found_instability
end

function test_stability(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, star_mass::Float64, sim_param::SimParam; ecc::AbstractVector{Float64}=zeros(length(P)))
    @assert length(P) == length(mass) == length(ecc)
    min_num_mutual_hill_radii = get_real(sim_param, "num_mutual_hill_radii")
    found_instability = false
    order = sortperm(P)
    a2 = semimajor_axis(P[order[1]], star_mass)
    for pl in 1:(length(P)-1)
        a1 = a2   # semimajor_axis(P[order[pl]],star_mass)
        a2 = semimajor_axis(P[order[pl+1]], star_mass)
        a = 0.5*(a1+a2)
        mu = (mass[order[pl]] + mass[order[pl+1]])/star_mass
        mutual_hill_radius = calc_hill_sphere(a, mu)
        e1 = ecc[order[pl]]
        e2 = ecc[order[pl+1]]
        if a2*(1-e2)-a1*(1+e1) < min_num_mutual_hill_radii*mutual_hill_radius
            found_instability = true
            break
        end
    end # loop over neighboring planet pairs within cluster
    return !found_instability
end

function is_period_ratio_near_resonance(period_ratio::Float64, sim_param::SimParam)
    resonance_width = get_real(sim_param, "resonance_width")
    resonance_width_factor = 1+resonance_width
    period_ratios_to_check = get_any(sim_param, "period_ratios_mmr", Array{Float64,1})
    result = false
    for period_ratio_mmr in period_ratios_to_check
        if period_ratio_mmr <= period_ratio <= period_ratio_mmr*resonance_width_factor
            result = true
            break
        end
    end
    return result
end

function calc_if_near_resonance(P::AbstractVector{Float64}, sim_param::SimParam)
    @assert issorted(P)   # TODO: OPT: Could remove once know it is safe
    result = falses(length(P))
    if length(P) >= 2
        for i in 1:(length(P)-1)
            if is_period_ratio_near_resonance(P[i+1]/P[i], sim_param)
                result[i] = true
                result[i+1] = true
            end # near mmr
        end # planets
    end # at least two planets
    return result
end


# Code to generate simple planetary systems.  For more sophisticated algorithms, see clustered model in He et al. 2019

function generate_planet_mass_from_radius_powerlaw(r::Float64, sim_param::SimParam)
  mr_power_index::Float64 = get_real(sim_param,"mr_power_index")
  mr_const::Float64 = get_real(sim_param,"mr_const")
  mr_max_mass::Float64 = get_real(sim_param,"mr_max_mass")
  m = mr_const*earth_mass*(r/earth_radius)^mr_power_index
  if m > mr_max_mass
     m = mr_max_mass
  end
  return m
end

function generate_planet_mass_from_radius_powerlaw(r::Float64, s::Star, o::Orbit, sim_param::SimParam)  # TODO USER SCI: This is importnat if you are using a stability criteria.  In that case, this should be replacde w/ better M-R relationship.  See Matthias's clustered  model for example.
  generate_planet_mass_from_radius_powerlaw(r,sim_param)
end

function generate_num_planets_poisson(lambda::Real, max_planets::Integer; min_planets::Integer = 0)
  ##### Note: this function produces odd behaviour if lambda < min_planets due to bugs in Distributions.Truncated(); for example, Distributions.Truncated(Distributions.Poisson(lambda),min_planets,max_planets) returns only values >=2 if min_planets=1 and lambda<~0.95
  if lambda < min_planets*1e-3
      return min_planets
  end
  bug_fixed = false # TODO OPT: true case should work, but Danley found bug in Distributions package.  Revert once fixed for speed.
  local n
  if bug_fixed
     d = Distributions.Truncated(Distributions.Poisson(lambda),min_planets,max_planets)
     n = rand(d)
  else
     if min_planets == 0
        min_planets = -1
     end
     d = Distributions.Truncated(Distributions.Poisson(lambda),min_planets,max_planets)
     n = rand(d)
     #=
     n = -1
     while !(min_planets<=n<=max_planets)
        n = rand(Distributions.Poisson(lambda))
     end
     =#
  end
  return n
end

function draw_truncated_poisson(lambda::Real; min::Integer=0, max::Integer=20, n::Integer=1)
    pmf = [(lambda^k*exp(-lambda))/factorial(k) for k in min:max]
    pmf ./= sum(pmf)
    cmf = zeros(Float64,length(pmf))
    for i in 1:length(pmf)
        for j in i:-1:1
            cmf[i] += pmf[j]
        end
    end
    result = Array{Int64}(undef,n)
    for i in 1:n
        u = rand()
        result[i] = findfirst(x-> x>=u,cmf) + min -1
    end
    return result
end

function generate_num_planets_poisson(s::Star, sim_param::SimParam)
  lambda::Float64 = exp(get_real(sim_param,"log_eta_pl"))
  max_tranets_in_sys::Int64 = get_int(sim_param,"max_tranets_in_sys")
  generate_num_planets_poisson(lambda,max_tranets_in_sys)
end

function generate_period_and_sizes_log_normal(s::Star, sim_param::SimParam; num_pl::Integer = 1)  # TODO USER SCI:  User should make sure planetary properties are being drawn appropriately for their scientific purposes
    mu_log_r::Float64 = get_real(sim_param,"mean_log_planet_radius")
    sigma_log_r::Float64 = get_real(sim_param,"sigma_log_planet_radius")
    mu_log_P::Float64 = get_real(sim_param,"mean_log_planet_period")
    sigma_log_P::Float64 = get_real(sim_param,"sigma_log_planet_period")
    min_period::Float64 = get_real(sim_param,"min_period")
    max_period::Float64 = get_real(sim_param,"max_period")
    min_radius::Float64 = get_real(sim_param,"min_radius")
    max_radius::Float64 = get_real(sim_param,"max_radius")
    max_draws::Int64 = 100

    if   sigma_log_r <= 0. || sigma_log_P<=0.
     println("# mu_log_r= ", mu_log_r, " sigma_log_r= ", sigma_log_r, " mu_log_P= ", mu_log_P, " sigma_log_P= ", sigma_log_P)
    end
    rdist = LogNormal(mu_log_r,sigma_log_r)
    Pdist = LogNormal(mu_log_P,sigma_log_P)
    #Rlist = rand(rdist,num_pl)
    #Plist = rand(Pdist,num_pl)
    #idx_keep = find(i->(min_radius<=Rlist[i]<=max_radius) && (min_period<=Plist[i]<=max_period), 1:num_pl )
    #return Plist[idx_keep], Rlist[idx_keep]  # replaced because want to return exactly num_pl.  Could use Truncated to restore above code.
    Rlist = zeros(num_pl)
    Plist = zeros(num_pl)
    for i in 1:num_pl
      j = 0
      while ! (min_radius<Rlist[i]<=max_radius) && j<max_draws
            Rlist[i] = rand(rdist)
            j+=1
      end
      if j>=max_draws
         println("# Struggled to draw size for: ",mu_log_r, " ", sigma_log_r)
      end
      j = 0
      while ! (min_period<Plist[i]<=max_period) && j<max_draws
            Plist[i] = rand(Pdist)
            j+=1
      end
      if j>=max_draws
         println("# Struggled to draw period for: ",mu_log_P, " ", sigma_log_P)
      end
    end
    return Plist, Rlist
end

function draw_power_law(n::Real, x0::Real, x1::Real, num_pl::Integer)
    if n != -1
        return ((x1^(n+1) - x0^(n+1)).*rand(num_pl) .+ x0^(n+1)).^(1/(n+1))
    else #if n == -1
        return x0*(x1/x0).^rand(num_pl)
    end
end

function draw_power_law!(out::Array, n::Real, x0::Real, x1::Real, num_pl::Integer)
    @assert length(out)==num_pl
    if n != -1
        return out .= ((x1^(n+1) - x0^(n+1)).*rand(num_pl) .+ x0^(n+1)).^(1/(n+1))
    else #if n == -1
        return out .= x0*(x1/x0).^rand(num_pl)
    end
end

function draw_broken_power_law(n1::Real, n2::Real, x0::Real, x1::Real, xb::Real, num_pl::Integer)
    #x0 and x1 are the lower and upper truncation limits, and xb is the break point, i.e. x0 <= xb <= x1 (all must be positive)
    #n1 and n2 are the power law indices between x0 and xb, and xb and x1, respectively (can be positive or negative)
    @assert(x0 <= xb <= x1)
    @assert(num_pl >= 1)

    u_draws = rand(num_pl) #'num_pl' draws from the uniform distribution between 0 and 1
    x_draws = zeros(num_pl)

    if (n1 != -1) & (n2 != -1)
        C1 = 1.0/(((xb^(n1+1) - x0^(n1+1))/(n1+1)) + ((xb^(n1-n2)*(x1^(n2+1) - xb^(n2+1)))/(n2+1))) #normalization constant
        ub = (C1*(xb^(n1+1) - x0^(n1+1)))/(n1+1) #break point in u, between 0 and 1
        for (i,u) in enumerate(u_draws)
            if u <= ub
                x_draws[i] = (((n1+1)*u)/C1 + x0^(n1+1))^(1/(n1+1))
            else #if u > ub
                x_draws[i] = (((n2+1)/(C1*xb^(n1-n2)))*(u - (C1*(xb^(n1+1) - x0^(n1+1)))/(n1+1)) + xb^(n2+1))^(1/(n2+1))
            end
        end
    elseif (n1 == -1) & (n2 != -1)
        C1 = 1.0/(log(xb/x0) + ((xb^(-1-n2))*(x1^(n2+1)) - 1)/(n2+1)) #normalization constant
        ub = C1*log(xb/x0) #break point in u, between 0 and 1
        for (i,u) in enumerate(u_draws)
            if u <= ub
                x_draws[i] = x0*exp(u/C1)
            else #if u > ub
                x_draws[i] = (((n2+1)/(C1*xb^(-1-n2)))*(u - C1*log(xb/x0)) + xb^(n2+1))^(1/(n2+1))
            end
        end
    elseif (n1 != -1) & (n2 == -1)
        C1 = 1.0/(((xb^(n1+1) - x0^(n1+1))/(n1+1)) + (xb^(n1+1))*log(x1/xb)) #normalization constant
        ub = (C1*(xb^(n1+1) - x0^(n1+1)))/(n1+1) #break point in u, between 0 and 1
        for (i,u) in enumerate(u_draws)
            if u <= ub
                x_draws[i] = (((n1+1)*u)/C1 + x0^(n1+1))^(1/(n1+1))
            else #if u > ub
                x_draws[i] = xb*exp((1/(C1*xb^(n1+1)))*(u - (C1*(xb^(n1+1) - x0^(n1+1)))/(n1+1)))
            end
        end
    else #if n1 == -1 and n2 == -1 (i.e. it is a single power-law with index of -1)
        for (i,u) in enumerate(u_draws)
            x_draws[i] = x0*exp(u*log(x1/x0))
        end
    end

    return x_draws
end

function generate_periods_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1)
    power_law_P::Float64 = get_real(sim_param,"power_law_P")
    min_period::Float64 = get_real(sim_param,"min_period")
    max_period::Float64 = get_real(sim_param,"max_period")
    Plist = draw_power_law(power_law_P,min_period,max_period, num_pl)
    return Plist
end

function generate_sizes_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1)
    power_law_r::Float64 = get_real(sim_param,"power_law_r")
    min_radius::Float64 = get_real(sim_param,"min_radius")
    max_radius::Float64 = get_real(sim_param,"max_radius")
    Rlist = draw_power_law(power_law_r,min_radius,max_radius, num_pl)
    return Rlist
end

function generate_sizes_broken_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1)
    power_law_r1::Float64 = get_real(sim_param,"power_law_r1")
    power_law_r2::Float64 = get_real(sim_param,"power_law_r2")
    min_radius::Float64 = get_real(sim_param,"min_radius")
    max_radius::Float64 = get_real(sim_param,"max_radius")
    break_radius::Float64 = get_real(sim_param,"break_radius")
    Rlist = draw_broken_power_law(power_law_r1,power_law_r2,min_radius,max_radius,break_radius, num_pl)
    return Rlist
end

function generate_period_and_sizes_power_law(s::Star, sim_param::SimParam; num_pl::Integer = 1)
    return (generate_periods_power_law(s, sim_param, num_pl=num_pl), generate_sizes_power_law(s, sim_param, num_pl=num_pl))
end





using SpecialFunctions

cdf_lognormal(x::Float64; μ::Float64=0., σ::Float64=1.) = 0.5*erfc(-(log(x) - μ)/(σ*sqrt(2.)))

function invert_cdf_lognormal(y::Float64; μ::Float64=0., σ::Float64=1.)
    @assert 0<=y<=1
    return exp(μ - sqrt(2.)*σ*erfcinv(2*y))
end

function cdf_power_law(x::Float64; x0::Float64, x1::Float64, α::Float64)
    @assert x0 <= x <= x1
    if α != -1
        return (x^(α+1) - x0^(α+1))/(x1^(α+1) - x0^(α+1))
    else # if α == -1
        return log(x/x0)/log(x1/x0)
    end
end

function invert_cdf_power_law(y::Float64; x0::Float64, x1::Float64, α::Float64)
    @assert x0 < x1
    if α != -1
        return (x0^(α+1) + (x1^(α+1) - x0^(α+1))*y)^(1/(α+1))
    else # if α == -1
        return x0*(x1/x0)^y
    end
end

function draw_segmented_uniform(segments::Array{Tuple{Float64,Float64},1})
    # Note: "segments" must be non-overlapping
    n_segs = length(segments) # number of segments
    @assert n_segs > 0
    for (i,seg) in enumerate(segments)
        @assert 0 <= seg[1] < seg[2] <= 1 # check that all segments are valid
        if i != n_segs
            @assert segments[i][2] < segments[i+1][1] # check to make sure there are no overlapping segments
        end
    end

    seg_lengths = [seg[2]-seg[1] for seg in segments]
    sum_lengths = sum(seg_lengths)

    u_draw = rand()
    for (i,seg) in enumerate(segments)
        # Map the union of the segments uniformly to the domain [0,1]:
        sum_l = sum(seg_lengths[1:i-1])
        if sum_l <= u_draw*sum_lengths < sum_l + seg_lengths[i]
            return rand(Uniform(seg[1], seg[2]))
        end
    end
end

function compute_unstable_regions_periods_given_planets(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, insert_pl_mass::Float64, star_mass::Float64, sim_param::SimParam; ecc::AbstractVector{Float64}=zeros(length(P)), insert_pl_ecc::Float64=0., use_mutualHill::Bool=true, verbose::Bool=false)
    @assert length(P) == length(mass) == length(ecc)
    min_num_mutual_hill_radii = get_real(sim_param, "num_mutual_hill_radii")
    order = sortperm(P)

    P_segments_unstable = Tuple{Float64,Float64}[]
    for pl in 1:length(P)
        a = semimajor_axis(P[order[pl]], star_mass)
        if use_mutualHill
            mu = (mass[order[pl]] + insert_pl_mass)/star_mass
            term = (min_num_mutual_hill_radii/2.)*(mu/3.)^(1//3)
            @assert 0 < 1. - insert_pl_ecc - term # must be true or else inequality flips sign
            a_lower, a_upper = a*(1. - ecc[order[pl]] - term)/(1. + insert_pl_ecc + term), a*(1. + ecc[order[pl]] + term)/(1. - insert_pl_ecc - term)
        else # use Hill radii
            mu = mass[order[pl]]/star_mass
            hill_radius = calc_hill_sphere(a, mu)
            a_lower, a_upper = a*(1. - ecc[order[pl]]) - hill_radius*min_num_mutual_hill_radii, a*(1. + ecc[order[pl]]) + hill_radius*min_num_mutual_hill_radii
        end
        a_lower = max(a_lower, 0.) # set a_lower=0 if negative
        @assert 0 <= a_lower < a_upper
        P_lower, P_upper = period_given_semimajor_axis(a_lower, insert_pl_mass+star_mass), period_given_semimajor_axis(a_upper, insert_pl_mass+star_mass)
        if verbose
            println("P_blocked: ", (P_lower, P_upper))
        end
        push!(P_segments_unstable, (P_lower, P_upper))
    end
    return P_segments_unstable
end

function compute_unstable_mutualHill_regions_periodscales_given_clusters(pl_per_cl::Vector{Int64}, P_cl::AbstractVector{Float64}, mass_cl::AbstractVector{Float64}, insert_cl_ρ::AbstractVector{Float64}, insert_cl_mass::AbstractVector{Float64}, star_mass::Float64, sim_param::SimParam; ecc_cl::AbstractVector{Float64}=zeros(length(P_cl)), insert_cl_ecc::AbstractVector{Float64}=zeros(length(insert_cl_mass)), verbose::Bool=false)
    n_cl = length(pl_per_cl)
    @assert sum(pl_per_cl) == length(P_cl) == length(mass_cl) == length(ecc_cl)
    @assert length(insert_cl_ρ) == length(insert_cl_mass) == length(insert_cl_ecc)
    min_num_mutual_hill_radii = get_real(sim_param, "num_mutual_hill_radii")

    # Get the inner-most and outer-most planets in the cluster we are trying to insert:
    order_insert_cl = sortperm(insert_cl_ρ) # ρ denotes unscaled periods in this function
    insert_pl_in_ρ, insert_pl_in_mass, insert_pl_in_ecc = insert_cl_ρ[order_insert_cl[1]], insert_cl_mass[order_insert_cl[1]], insert_cl_ecc[order_insert_cl[1]]
    insert_pl_out_ρ, insert_pl_out_mass, insert_pl_out_ecc = insert_cl_ρ[order_insert_cl[end]], insert_cl_mass[order_insert_cl[end]], insert_cl_ecc[order_insert_cl[end]]
    max_ratio = insert_pl_out_ρ/insert_pl_in_ρ

    Pc_segments_unstable = Tuple{Float64,Float64}[]
    pl_start = 1
    pl_stop = 0
    for c in 1:n_cl
        pl_stop += pl_per_cl[c]
        order = (pl_start - 1) .+ sortperm(P_cl[pl_start:pl_stop])

        # Compare the inner-most planet in this cluster to the outer-most planet in the cluster we are trying to fit in to calculate the lower bound for the unstable region:
        a_in = semimajor_axis(P_cl[order[1]], star_mass)
        mu = (mass_cl[order[1]] + insert_pl_out_mass)/star_mass
        term = (min_num_mutual_hill_radii/2.)*(mu/3.)^(1//3)
        a_lower = a_in*(1. - ecc_cl[order[1]] - term)/(1. + insert_pl_out_ecc + term)
        a_lower = max(a_lower, 0.) # set a_lower=0 if negative
        P_lower = period_given_semimajor_axis(a_lower, insert_pl_out_mass+star_mass)

        # Compare the outer-most planet in this cluster to the inner-most planet in the cluster we are trying to fit in to calculate the upper bound for the unstable region:
        a_out = semimajor_axis(P_cl[order[end]], star_mass)
        mu = (mass_cl[order[end]] + insert_pl_in_mass)/star_mass
        term = (min_num_mutual_hill_radii/2.)*(mu/3.)^(1//3)
        @assert 0 < 1. - insert_pl_in_ecc - term # must be true or else inequality flips sign
        a_upper = a_out*(1. + ecc_cl[order[end]] + term)/(1. - insert_pl_in_ecc - term)
        P_upper = period_given_semimajor_axis(a_upper, insert_pl_out_mass+star_mass)

        # Now, convert P_lower and P_upper to bounds for the period scale (which depend on the unscaled periods of the cluster we are trying to insert):
        Pc_lower, Pc_upper = P_lower/insert_pl_out_ρ, P_upper/insert_pl_in_ρ # Pc_lower should decrease if insert_pl_out_ρ > 1, and likewise Pc_upper should increase if insert_pl_in_ρ < 1
        if verbose
            println("Pc_blocked: ", (Pc_lower, Pc_upper))
        end
        push!(Pc_segments_unstable, (Pc_lower, Pc_upper))

        pl_start += pl_per_cl[c]
    end
    Pc_segments_unstable = sort(Pc_segments_unstable)

    Pc_segments_unstable_nonoverlapping = Tuple{Float64,Float64}[]
    x_start, x_stop = Pc_segments_unstable[1]
    for i in 1:length(Pc_segments_unstable)-1
        if x_stop > Pc_segments_unstable[i+1][1] # segments overlap
            x_stop = Pc_segments_unstable[i+1][2]
        else # segments do not overlap
            if Pc_segments_unstable[i+1][1]/x_stop > max_ratio
                push!(Pc_segments_unstable_nonoverlapping, (x_start, x_stop))
                x_start, x_stop = Pc_segments_unstable[i+1]
            else # segments do not overlap, but there is not enough room between the two segments to insert the cluster
                x_stop = Pc_segments_unstable[i+1][2]
            end
        end
    end
    push!(Pc_segments_unstable_nonoverlapping, (x_start, x_stop))
    if verbose
        println("Pc_blocked with no room: ", Pc_segments_unstable_nonoverlapping)
    end
    return Pc_segments_unstable_nonoverlapping
end

function compute_blocked_regions_non_overlapping(segments_blocked::Array{Tuple{Float64,Float64},1})
    # NOTE: this function assumes that the input segments are sorted in increasing order!
    # Make a list of non-overlapping segments:
    segments_blocked_nonoverlapping = Tuple{Float64,Float64}[]
    x_start, x_stop = segments_blocked[1]
    for i in 1:length(segments_blocked)-1
        if x_stop > segments_blocked[i+1][1] # segments overlap
            x_stop = segments_blocked[i+1][2]
        else # segments do not overlap
            push!(segments_blocked_nonoverlapping, (x_start, x_stop))
            x_start, x_stop = segments_blocked[i+1]
        end
    end
    push!(segments_blocked_nonoverlapping, (x_start, x_stop))
    return segments_blocked_nonoverlapping
end

function compute_allowed_regions_cdf_lognormal(segments_blocked::Array{Tuple{Float64,Float64},1}; μ::Float64=0., σ::Float64=1., x_min::Float64=0., x_max::Float64=Inf)
    @assert 0 <= x_min < x_max <= Inf

    # Make a list of non-overlapping segments in period deemed unstable:
    segments_blocked_nonoverlapping = compute_blocked_regions_non_overlapping(segments_blocked)

    # Make a list of non-overlapping, allowed segments in the cdf:
    cdf_segments_allowed = Tuple{Float64,Float64}[]
    n_seg_max = length(segments_blocked_nonoverlapping)+1 # maximum number of allowed segments if all the unstable segments fit between x_min and x_max
    for i in 1:n_seg_max
        if i==1 && x_min < segments_blocked_nonoverlapping[i][1]
            cdf_start = cdf_lognormal(x_min; μ=μ, σ=σ)
            cdf_stop = cdf_lognormal(segments_blocked_nonoverlapping[i][1]; μ=μ, σ=σ)
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        elseif i==n_seg_max && x_max > segments_blocked_nonoverlapping[i-1][2]
            cdf_start = cdf_lognormal(segments_blocked_nonoverlapping[i-1][2]; μ=μ, σ=σ)
            cdf_stop = cdf_lognormal(x_max; μ=μ, σ=σ)
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        elseif 1 < i < n_seg_max
            cdf_start = cdf_lognormal(segments_blocked_nonoverlapping[i-1][2]; μ=μ, σ=σ)
            cdf_stop = cdf_lognormal(segments_blocked_nonoverlapping[i][1]; μ=μ, σ=σ)
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        end
    end

    return cdf_segments_allowed
end

function compute_allowed_regions_cdf_power_law(segments_blocked::Array{Tuple{Float64,Float64},1}; x0::Float64, x1::Float64, α::Float64)
    @assert x0 < x1

    # Make a list of non-overlapping segments in period scale deemed unstable:
    segments_blocked_nonoverlapping = compute_blocked_regions_non_overlapping(segments_blocked)

    # Make a list of non-overlapping, allowed segments in the cdf:
    cdf_segments_allowed = Tuple{Float64,Float64}[]
    n_seg_max = length(segments_blocked_nonoverlapping)+1 # maximum number of allowed segments if all the unstable segments fit between x0 and x1
    for i in 1:n_seg_max
        if i==1 && x0 < segments_blocked_nonoverlapping[i][1]
            cdf_start = 0.
            cdf_stop = cdf_power_law(segments_blocked_nonoverlapping[i][1]; x0=x0, x1=x1, α=α)
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        elseif i==n_seg_max && x1 > segments_blocked_nonoverlapping[i-1][2]
            cdf_start = cdf_power_law(segments_blocked_nonoverlapping[i-1][2]; x0=x0, x1=x1, α=α)
            cdf_stop = 1.
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        elseif 1 < i < n_seg_max
            cdf_start = cdf_power_law(segments_blocked_nonoverlapping[i-1][2]; x0=x0, x1=x1, α=α)
            cdf_stop = cdf_power_law(segments_blocked_nonoverlapping[i][1]; x0=x0, x1=x1, α=α)
            push!(cdf_segments_allowed, (cdf_start, cdf_stop))
        end
    end

    return cdf_segments_allowed
end

function draw_lognormal_allowed_regions(segments_blocked::Array{Tuple{Float64,Float64},1}; μ::Float64=0., σ::Float64=1., x_min::Float64=0., x_max::Float64=Inf, verbose=false)
    cdf_segments_allowed = compute_allowed_regions_cdf_lognormal(segments_blocked; μ=μ, σ=σ, x_min=x_min, x_max=x_max)
    if length(cdf_segments_allowed) == 0
        if verbose
            println("No allowed regions left to draw from; returning NaN.")
        end
        return NaN
    else
        y = draw_segmented_uniform(cdf_segments_allowed)
        return invert_cdf_lognormal(y; μ=μ, σ=σ)
    end
end

function draw_power_law_allowed_regions(segments_blocked::Array{Tuple{Float64,Float64},1}; x0::Float64, x1::Float64, α::Float64, verbose=false)
    cdf_segments_allowed = compute_allowed_regions_cdf_power_law(segments_blocked; x0=x0, x1=x1, α=α)
    if length(cdf_segments_allowed) == 0
        if verbose
            println("No allowed regions left to draw from; returning NaN.")
        end
        return NaN
    else
        y = draw_segmented_uniform(cdf_segments_allowed)
        return invert_cdf_power_law(y; x0=x0, x1=x1, α=α)
    end
end

function draw_period_lognormal_allowed_regions(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, insert_pl_mass::Float64, star_mass::Float64, sim_param::SimParam; μ::Float64=0., σ::Float64=1., x_min::Float64=0., x_max::Float64=Inf, ecc::AbstractVector{Float64}=zeros(length(P)), insert_pl_ecc::Float64=0., use_mutualHill::Bool=true, verbose::Bool=false)
    @assert length(P) == length(mass) == length(ecc)
    @assert 0 <= insert_pl_mass < star_mass
    @assert all(P .> 0)
    @assert all(mass .>= 0)
    @assert all(1 .> ecc .>= 0)
    @assert 0 <= insert_pl_ecc < 1

    if length(P) > 0
        P_segments_unstable = compute_unstable_regions_periods_given_planets(P, mass, insert_pl_mass, star_mass, sim_param; ecc=ecc, insert_pl_ecc=insert_pl_ecc, use_mutualHill=use_mutualHill, verbose=verbose)
        P_draw = draw_lognormal_allowed_regions(P_segments_unstable; μ=μ, σ=σ, x_min=x_min, x_max=x_max, verbose=verbose)
    else
        P_draw = invert_cdf_lognormal(rand(Uniform(cdf_lognormal(x_min; μ=μ, σ=σ), cdf_lognormal(x_max; μ=μ, σ=σ))); μ=μ, σ=σ)
    end
    return P_draw
end

draw_period_lognormal_allowed_regions_Hill(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, insert_pl_mass::Float64, star_mass::Float64, sim_param::SimParam; μ::Float64=0., σ::Float64=1., x_min::Float64=0., x_max::Float64=Inf, ecc::AbstractVector{Float64}=zeros(length(P)), insert_pl_ecc::Float64=0., verbose::Bool=false) = draw_period_lognormal_allowed_regions(P, mass, insert_pl_mass, star_mass, sim_param; μ=μ, σ=σ, x_min=x_min, x_max=x_max, ecc=ecc, insert_pl_ecc=insert_pl_ecc, use_mutualHill=false, verbose=verbose)
draw_period_lognormal_allowed_regions_mutualHill(P::AbstractVector{Float64}, mass::AbstractVector{Float64}, insert_pl_mass::Float64, star_mass::Float64, sim_param::SimParam; μ::Float64=0., σ::Float64=1., x_min::Float64=0., x_max::Float64=Inf, ecc::AbstractVector{Float64}=zeros(length(P)), insert_pl_ecc::Float64=0., verbose::Bool=false) = draw_period_lognormal_allowed_regions(P, mass, insert_pl_mass, star_mass, sim_param; μ=μ, σ=σ, x_min=x_min, x_max=x_max, ecc=ecc, insert_pl_ecc=insert_pl_ecc, use_mutualHill=true, verbose=verbose)

function draw_periodscale_power_law_allowed_regions_mutualHill(pl_per_cl::Vector{Int64}, P_cl::AbstractVector{Float64}, mass_cl::AbstractVector{Float64}, insert_cl_ρ::AbstractVector{Float64}, insert_cl_mass::AbstractVector{Float64}, star_mass::Float64, sim_param::SimParam; x0::Float64, x1::Float64, α::Float64, ecc_cl::AbstractVector{Float64}=zeros(length(P_cl)), insert_cl_ecc::AbstractVector{Float64}=zeros(length(insert_cl_mass)), verbose::Bool=false)
    @assert x0 < x1
    @assert sum(pl_per_cl) == length(P_cl) == length(mass_cl) == length(ecc_cl)
    @assert all(P_cl .> 0)
    @assert all(mass_cl .> 0)
    @assert all(0 .<= ecc_cl .< 1)
    @assert length(insert_cl_ρ) == length(insert_cl_mass) == length(insert_cl_ecc)
    @assert all(0 .< insert_cl_ρ)
    @assert all(0 .< insert_cl_mass .< star_mass)
    @assert all(0 .<= insert_cl_ecc .< 1)

    if length(P_cl) > 0
        Pc_segments_unstable = compute_unstable_mutualHill_regions_periodscales_given_clusters(pl_per_cl, P_cl, mass_cl, insert_cl_ρ, insert_cl_mass, star_mass, sim_param; ecc_cl=ecc_cl, insert_cl_ecc=insert_cl_ecc, verbose=verbose)
        Pc_draw = draw_power_law_allowed_regions(Pc_segments_unstable; x0=x0, x1=x1, α=α, verbose=verbose)
    else
        Pc_draw = invert_cdf_power_law(rand(); x0=x0, x1=x1, α=α)
    end
    return Pc_draw
end





function TruncatedUpper(d::Distributions.UnivariateDistribution, u::Float64)
    zero(u) < u || error("lower bound should be less than upper bound.")
    lcdf = zero(u)
    ucdf = isinf(u) ? one(u) : cdf(d, u)
    tp = ucdf - lcdf
    Distributions.Truncated{typeof(d),Distributions.value_support(typeof(d))}(d, zero(u), u, lcdf, ucdf, tp, log(tp))
end


function generate_e_omega_rayleigh_direct(sigma_hk::Float64; max_e::Float64 = 1.0)
  @assert(0<max_e<=1.0)
  ecc::Float64 = rand( TruncatedUpper(Rayleigh(sigma_hk),max_e) )
  w::Float64 = 2pi*rand()
  return ecc, w
end

function generate_e_omega_rayleigh_two_gaussians(sigma_hk::Float64; max_e::Float64 = 1.0)
  @assert(0<max_e<=1.0)
  h = k = 1.0
  while h*h+k*k >= max_e*max_e
    h = sigma_hk*randn()
    k = sigma_hk*randn()
  end
  ecc::Float64 = sqrt(h*h+k*k)
  #w::Float64 = atan2(k,h)
  w::Float64 = atan(k,h)
  return ecc, w
end

function generate_e_omega_rayleigh(sigma_hk::Float64; max_e::Float64 = 1.0)
   if max_e > sigma_hk
      return generate_e_omega_rayleigh_two_gaussians(sigma_hk,max_e=max_e)
   else
      return generate_e_omega_rayleigh_direct(sigma_hk,max_e=max_e)
   end
end

function generate_e_omega_rayleigh(sim_param::SimParam; max_e::Float64 = 1.0)
  sigma_hk::Float64 = get_real(sim_param,"sigma_hk")
  generate_e_omega_rayleigh(sigma_hk, max_e=max_e)
end

function map_square_to_triangle(r1::Float64, r2::Float64, A::Vector{Float64}, B::Vector{Float64}, C::Vector{Float64})
    #This function takes in a point (r1,r2) in the unit square (i.e. r1,r2 in [0,1]) and maps it to a point P=(x,y) in the triangle defined by vertices A,B,C
    #If r1,r2 are uniformly drawn in [0,1], then the point P=(x,y) is also uniformly drawn in the triangle; see http://www.cs.princeton.edu/~funk/tog02.pdf (Section 4.2) for a reference

    @assert 0. <= r1 <= 1.
    @assert 0. <= r2 <= 1.
    P = (1. - sqrt(r1)) .* A + (sqrt(r1)*(1. - r2)) .* B + (sqrt(r1)*r2) .* C
    return P
end

function generate_planetary_system_hardcoded_example(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # in this version we specify fixed functions that are known at compile time, allowing for additional optimizations (~0.6 second faster per Kepler catalog out of ~3.6 sec on my laptop w/ 1 core)
  generate_planet_mass_from_radius = generate_planet_mass_from_radius_powerlaw
  generate_num_planets = generate_num_planets_poisson
  generate_period_and_sizes = generate_period_and_sizes_log_normal

  generate_e_omega =  generate_e_omega_rayleigh

  # generate_star = get_function(sim_param,"generate_star")
  # star::StarAbstract = generate_star(sim_param)
  num_pl::Int64 = generate_num_planets(star, sim_param)

  if( num_pl==0 )
    return PlanetarySystem(star)
  else
    (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
    idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred

    min_a_in_rstar = 2.0
    min_P_orbit = day_in_year*sqrt((min_a_in_rstar*star.radius*rsol_in_au)^3 / star.mass) # minimum semi-major axis of two stellar radii
    idx = idx[Plist[idx] .> min_P_orbit]
    if( length(idx)==0 )
        return PlanetarySystem(star)
    end

    pl = Array{Planet}(length(idx))
    orbit = Array{Orbit}(length(idx))
    a = map(i->semimajor_axis(Plist[i],star.mass),idx)
    max_e = ones(length(idx))
    max_e_factor = 0.999 # A factor just less than 1 to prevent numerical issues with near-crossing orbits
    if length(a)>=2
       for i in 1:length(a)
         if i==1
           max_e[i] = max_e_factor*(1-a[i]/a[i+1])/(1+a[i]/a[i+1])
         elseif i==length(a)
           max_e[i] = max_e_factor*(1-a[i-1]/a[i])/(1+a[i-1]/a[i])
         else
           max_e[i] = max_e_factor*min( (1-a[i]/a[i+1])/(1+a[i]/a[i+1]), (1-a[i-1]/a[i])/(1+a[i-1]/a[i]) )
         end
       end
    end
    for i in 1:length(idx)
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      Rpl = Rlist[idx[i]]
      (ecc::Float64,  omega::Float64) = generate_e_omega(sim_param, max_e=max_e[i])
      incl::Float64 = acos(rand())
      orbit[i] = Orbit(P,ecc,incl,omega,2pi*rand(),2pi*rand())
      mass::Float64 = generate_planet_mass_from_radius(Rpl, star, orbit[i], sim_param)
      pl[i] = Planet( Rpl,  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

function generate_planetary_system_uncorrelated_incl(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # load functions to use for drawing parameters
   generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")
   generate_num_planets = get_function(sim_param,"generate_num_planets")
   generate_period_and_sizes = get_function(sim_param,"generate_period_and_sizes")
   generate_e_omega = get_function(sim_param,"generate_e_omega")

  #  generate_star = get_function(sim_param,"generate_star")
  # star::StarAbstract = generate_star(sim_param)
  num_pl::Int64 = generate_num_planets(star, sim_param)::Int64
  sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0

  if( num_pl==0 )
    return PlanetarySystem(star)::PlanetarySystem
  else
     (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
    idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred

    min_a_in_rstar = 2.0
    min_P_orbit = day_in_year*sqrt((min_a_in_rstar*star.radius*rsol_in_au)^3 / star.mass) # minimum semi-major axis of two stellar radii
    idx = idx[Plist[idx] .> min_P_orbit]
    if( length(idx)==0 )
        return PlanetarySystem(star)
    end

    pl = Array{Planet}(undef,length(idx))
    orbit = Array{Orbit}(undef,length(idx))
    a = map(i->semimajor_axis(Plist[i],star.mass),idx)
    max_e = ones(length(idx))
    max_e_factor = 0.999 # A factor just less than 1 to prevent numerical issues with near-crossing orbits
    if length(a)>=2
       for i in 1:length(a)
          if i==1
            max_e[i] = max_e_factor*(1-a[i]/a[i+1])/(1+a[i]/a[i+1])
          elseif i==length(a)
            max_e[i] = max_e_factor*(1-a[i-1]/a[i])/(1+a[i-1]/a[i])
          else
            max_e[i] = max_e_factor*min( (1-a[i]/a[i+1])/(1+a[i]/a[i+1]), (1-a[i-1]/a[i])/(1+a[i-1]/a[i]) )
          end
        end
    end
    for i in 1:length(idx)
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      Rpl = Rlist[idx[i]]
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc::Float64,  omega::Float64) = generate_e_omega(sim_param, max_e=max_e[i])
      incl::Float64 = acos(rand())
      orbit[i] = Orbit(P,ecc,incl,omega,2pi*rand(),2pi*rand())
      # set!(orbit[idx[i]],P,ecc,incl,omega,2pi*rand(),2pi*rand())
      mass::Float64 = generate_planet_mass_from_radius(Rpl, star, orbit[i], sim_param)
      pl[i] = Planet( Rpl,  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

# This version generates more systems roughly near a common plane, but until incorporate CORBITS data, ABC can't match input param
function generate_planetary_system_simple(star::StarAbstract, sim_param::SimParam; verbose::Bool = false)
  # load functions to use for drawing parameters
    generate_planet_mass_from_radius = get_function(sim_param,"generate_planet_mass_from_radius")
    generate_num_planets = get_function(sim_param,"generate_num_planets")
  #  generate_num_planets = generate_num_planets_christiansen
  #  generate_period_and_sizes = generate_period_and_sizes_christiansen
    generate_period_and_sizes = get_function(sim_param,"generate_period_and_sizes")
    generate_e_omega = get_function(sim_param,"generate_e_omega")
    sigma_incl = deg2rad(get_real(sim_param,"sigma_incl"))

  #   generate_star = get_function(sim_param,"generate_star")
  #  star::StarAbstract = generate_star(sim_param)
   num_pl = generate_num_planets(star, sim_param)::Int64
   sigma_ecc::Float64 = haskey(sim_param,"sigma_hk") ? get_real(sim_param,"sigma_hk") : 0.0

  if( num_pl==0 )
    return PlanetarySystem(star)
  else
    (Plist::Vector{Float64}, Rlist::Vector{Float64}) = generate_period_and_sizes(star, sim_param, num_pl=num_pl)
    idx = sortperm(Plist)                   # TODO OPT: Check to see if sorting is significant time sink.  If so, it might could be deferred
    incl_sys = acos(rand())

    min_a_in_rstar = 2.0
    min_P_orbit = day_in_year*sqrt((min_a_in_rstar*star.radius*rsol_in_au)^3 / star.mass) # minimum semi-major axis of two stellar radii
    idx = idx[Plist[idx] .> min_P_orbit]
    if( length(idx)==0 )
        return PlanetarySystem(star)
    end

    pl = Array{Planet}(undef,length(idx))
    orbit = Array{Orbit}(undef,length(idx))
    for i in 1:length(idx)
      # if verbose   println("i=",i," idx=",idx," Plist=",Plist[idx] );     end
      P = Plist[idx[i]]
      Rpl = Rlist[idx[i]]
      if haskey(sim_param,"sigma_hk_one") && haskey(sim_param,"sigma_hk_multi")
         sigma_ecc = num_pl == 1 ? get_real(sim_param,"sigma_hk_one") : get_real(sim_param,"sigma_hk_multi")
      end
      (ecc,  omega) = generate_e_omega(sim_param)::Tuple{Float64,Float64}
      incl_mut = sigma_incl*sqrt(randn()^2+randn()^2) # rand(Distributions.Rayleigh(sigma_incl)) # sigma_incl*randn()
      asc_node = 2pi*rand()
      mean_anom = 2pi*rand()
      #incl = incl_sys + sigma_incl*randn()
      incl =  incl_mut!=zero(incl_mut) ? acos( cos(incl_sys)*cos(incl_mut) + sin(incl_sys)*sin(incl_mut)*cos(asc_node) ) : incl_sys
      orbit[i] = Orbit(P,ecc,incl,omega,asc_node,mean_anom)
      mass = generate_planet_mass_from_radius(Rpl, star, orbit[i], sim_param)::Float64
      pl[i] = Planet( Rpl,  mass )
    end
  return PlanetarySystem(star,pl,orbit)
  end
end

function test_planetary_system_constructors(sim_param::SimParam)
  generate_star = get_function(sim_param,"generate_star")
  star = generate_star(sim_param)
  empty_sys = PlanetarySystem(star)
  earth = Planet(earth_radius,earth_mass)
  earth_orbit = Orbit(365.2425,0.0167,0.5*pi,0.0,0.0,0.0)
  solar_sys = PlanetarySystem(star, earth,earth_orbit)
  m = generate_planet_mass_from_radius_powerlaw(0.02,star,earth_orbit,sim_param)/earth_mass
  generate_planetary_system_simple(star,sim_param,verbose=true)
end
