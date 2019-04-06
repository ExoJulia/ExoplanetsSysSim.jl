using ExoplanetsSysSim

start_dir = pwd()
data_dir = joinpath(dirname(pathof(ExoplanetsSysSim)),"..","data")
pkg_dev_home = joinpath(dirname(pathof(ExoplanetsSysSim)),"..")

cd(pkg_dev_home)
if true
   if !isfile(".gitmodules")
      println("# Adding git@github.com:ExoJulia/SysSimData.git a submodule in the data directory...")
      flush(stdout)
      run(`git submodule add git@github.com:ExoJulia/SysSimData.git data`)
   end
   println("# Initializing the data submodule...")
   flush(stdout)
   run(`git submodule init`)
   run(`git submodule update`)
   cd("data")
   println("# Initializing data's submodules...")
   flush(stdout)
   run(`git submodule init`)
   println("# git lfs pull just in case binary files not downloaded already...")
   flush(stdout)
   run(`git lfs pull`)
   cd("..")
end
println("# Recursively updating submodules...")
flush(stdout)
run(`git submodule update --recursive`)

if data_dir != joinpath(pkg_dev_home,"data")
  cd(pkg_dev_home)
  symlink(data_dir,"data")
end

cd(start_dir)

