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
