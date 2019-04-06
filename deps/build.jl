using ExoplanetsSysSim

start_dir = pwd()
data_dir_home = joinpath(dirname(pathof(ExoplanetsSysSim)),"..")
pkg_dev_home = joinpath(dirname(pathof(ExoplanetsSysSim)),"..")

cd(pkg_dev_home)
#= Disabled to prevent merge conflicts
println("# Pulling most recent version, just to be sure.")
flush(stdout)
run(`git pull`)
=#

# Disabled, so can create modules manually in another directory and make data a symlink 
# include("init_modules.jl")

cd(start_dir)

