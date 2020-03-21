## ExoplanetsSysSim/src/orbit.jl
## (c) 2015 Eric B. Ford

#export Orbit

struct Orbit     
  P::Float64             # days             # QUERY:  Should we store P or a here?
  ecc::Float64
  incl::Float64          # radians; inclination relative to sky of observer
  omega::Float64         # radians
  asc_node::Float64      # radians
  mean_anom::Float64     # radians          # QUERY:  Should we store t_0 or mean_anom here?
end

#Orbit() = Orbit(0.0,0.0,0.0,0.0,0.0,0.0)   # Comment out, so don't accidentally have invalid orbits

# This will only work if Orbit were mutable.  Is that better or worse?  Let's test and see....
function set!(o::Orbit, P::Float64, e::Float64, i::Float64, w::Float64, asc_node::Float64, M::Float64)
 o.P = P
 o.ecc = e
 o.incl = i
 o.omega = w
 o.asc_node = asc_node
 o.mean_anom = M
 return o
end


function test_orbit_constructors()
  #orb = Orbit()
  orb = Orbit(1.0, 0.03, 0.5*pi, 0.0, 0.0, pi)
  if !isimmutable(orb)
     set!(orb,1.0, 0.03, 0.5*pi, 0.0, 0.0, pi)
  end
end




"""

Solve for `i_m` in the spherical law of Cosines: cos(i) = cos(i_ref)*cos(i_m) + sin(i_ref)*sin(i_m)*cos(Ω).
"""
#=
function calc_inclmut_to_refplane(o::Orbit, sp::SystemPlane)
    i = o.incl
    Ω = o.asc_node
    i_ref = sp.incl

    local i_m::Float64

    if cos(i_ref) ≈ -cos(i)
        i_m = π
    else
        cond1 = -cos(Ω)*sin(i_ref)*sqrt(cos(Ω)^2 * sin(i_ref)^2 - cos(i)^2 + cos(i_ref)^2) + cos(Ω)^2 * sin(i_ref)^2 + cos(i)*cos(i_ref) + cos(i_ref)^2
        cond2 = cos(Ω)*sin(i_ref)*sqrt(cos(Ω)^2 * sin(i_ref)^2 - cos(i)^2 + cos(i_ref)^2) + cos(Ω)^2 * sin(i_ref)^2 + cos(i)*cos(i_ref) + cos(i_ref)^2
        @info(cond1, cond2)
        if(cond1 != 0)
            numer = cos(Ω)*sin(i_ref) - sqrt( cos(Ω)^2 * sin(i_ref)^2 - cos(i)^2 + cos(i_ref)^2 )
        elseif(cond2 != 0)
            numer = cos(Ω)*sin(i_ref) + sqrt( cos(Ω)^2 * sin(i_ref)^2 - cos(i)^2 + cos(i_ref)^2 )
        else
            @warn("Unspecified case when computing mutual inclination!")
        end
        denom = cos(i) + cos(i_ref)
        i_m = 2*atan(numer/denom)
    end

    @assert(cos(i) ≈ cos(i_ref)*cos(i_m) + sin(i_ref)*sin(i_m)*cos(Ω))
    return i_m
end
=#
