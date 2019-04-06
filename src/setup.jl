## ExoplanetsSysSim/src/setup.jl
## (c) 2015 Eric B. Ford


if VERSION >=v"0.7.0-"
   using Pkg
   #using Libdl
end


# How to install ExoplanetsSysSim package
# Pkg.add(PackageSpec(url="git@github.com:eford/ExoplanetsSysSim.jl.git"))

#=  
# Is this still needed, now that we've moved to github?
# If so, this would need to be updated for Julia v0.7.0
# Since bitbucket messed up capitalization of package name
if ! isdir( joinpath(Pkg.devdir(),"ExoplanetsSysSim") )
     symlink( joinpath(Pkg.dir(),"exoplanetssyssim"), joinpath(Pkg.dir(),"ExoplanetsSysSim") )
end
=#

try
  Pkg.add(PackageSpec(url="git@github.com:eford/ABC.jl.git"))
catch
  warn("Attempted to install ABC.jl package, but was not successful.")
  warn("While most of SysSim will still work, some functionality will not be avaliable unless you install ABC correctly.")
end

try
  Pkg.add(PackageSpec(url="git@github.com:jbrakensiek/CORBITS.git"))

  # Compile CORBITS library and put it somewhere we can find
  cd(joinpath(Pkg.devdir(),"CORBITS"))
  run(`make lib`)
  cd(homedir())
  if !is_windows()
     symlink( joinpath(Pkg.devdir(),"CORBITS","libcorbits.so"), joinpath(Pkg.devdir(),"ExoplanetsSysSim","libcorbits.so") )
  else
     cp( joinpath(Pkg.devdir(),"CORBITS","libcorbits.so"), joinpath(Pkg.devdir(),"ExoplanetsSysSim","libcorbits.so") )
  end

catch
  warn("Attempted to install CORBITS.jl package, but was not successful.")
  warn("While most of SysSim will still work, some functionality will not be avaliable unless you install CORBITS correctly.")
end

