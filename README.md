# ExoplanetsSysSim
Welcome to the ExoplanetsSysSim package for generating planetary systems and simulating observations of those systems with a transit survey.  Currently, SysSim focuses on NASA's Kepler mission, but we've aimed to develop a generic framework that can be applied to other surveys (e.g., K2, TESS, PLATO, LSST, etc.).

# How to install ExoplanetsSysSim:
* Make sure you have julia (v0.7 or greater) installed
* Make sure you have git and [git-lfs](https://git-lfs.github.com/) installed.
If you're using ICS-ACI, then do can "instal" git-lfs by
```sh
export PATH=/gpfs/group/dml129/default/sw/git-lfs:$PATH
module load git
```
* If you want to use ssh keys instead of https authentication (to minimize typing your github password), then:
  * Setup a local ssh key using ssh-keygen
  * Tell Github about your ssh key:  Person Icon (upper right), Settings, SSH & GPG keys, New SSH Key.  Entry a name in the title box and paste the contents of `cat ~/.ssh/id_rsa.pub` into the "Key" box. Add SSH Key.
* Make sure that your JULIA_DEPOT_PATH (~/.julia by default) does not include an old version of CORBITS or ExopalnetsSysSim.  If this is your first time using julia v1.0, then you probably don't need to do anything.  Otherwise, I see two ways to do this:
   - One way to avoid conflicts is to delete the old repos.  But first, if there's any chance that you might have things in your current CORBITS or ExoplanetsSysSim repots that you want to keep, then make a backup copy of thsoe repos before proceeding.  Then 'Pkg.rm("CORBITS"); Pkg.rm("ExoplanetsSysSim"); Pkg.gc()' and rm -rf both from the dev subdirectory of your JULIA_DEPOT_PATH (~/.julia by default).  Warning:  Sometimes Julia manages to keep these around despite my efforts to delete them.
   - Another way to avoid conflicts with old versions is to sepcify a new JULIA_DEPOT_PATH.  If you go this route, then you'll need to make sure that this environment variable is set to the desired depot in each of your future shell sessions.
```sh
export JULIA_DEPOT_PATH=~/.julia_clean
```

* Run julia and install eford's version of the CORBITS.jl Julia package
```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/eford/CORBITS.jl"))
```
* Run julia and install the ExoplanetsSysSim repo as a Julia package using the following command.
```julia
Pkg.develop(PackageSpec(url="git@github.com:ExoJulia/ExoplanetsSysSim.jl.git"))
```
* Get the required data files to the data directory
   - Ideally, the deps/build.jl script will have automatically setup submodules, populating your data directory
   - If not, then you can manually create a data directory by cloning https://github.com/ExoJulia/SysSimData
* If creating the data directory manually  somewhere outside of the ExoplanetsSysSim repo, then make data files avaliable to SysSim, either by
   - Addin a symlink from whever you cloned SysSimData to '${JULIA_DEPOT_PATH}/dev/ExoplanetsSysSim'
   - Overriding the default file paths to point to whereever you placed the binary input files

* TODO CREATE TESTS.  Optionally, run some tests, e.g.
```julia
using ExoplanetsSysSim
include(joinpath(dirname(pathof(ExoplanetsSysSim)),"..","test","runtests.jl"))
```
* To use SysSim:
   - Create your own repository containing code that will call ExoplanetsSysSim
   - Make it a Julia project by adding dependancies, including CORBITS
   - Make your project depend on your development branch of ExoplanetsSysSim
```julia
using Pkg
Pkg.activate(".")
Pkg.add(PackageSpec(url="https://github.com/eford/CORBITS.jl"))
Pkg.develop(PackageSpec(url="git@github.com:ExoJulia/ExoplanetsSysSim.jl.git"))
```
   - Have your project code load ExoplanetsSysSim and use it
```julia
using ExoplanetsSysSim
```
* Write papers and cite relevant publications


# Team:
## Developers:
  * Eric Ford
  * Matthias He
  * Danley Hsu
  * Darin Ragozzine
## Other Contributors/Consultants:
  * Robert Morehead
  * Keir Ashby
  * Jessi Cisewski
  * Chad Schafer
  * Tom Loredo
  * Robert Wolpert

# Acknowledgements:
* NASA
  * Kepler Mission
  * Kepler Science Team
  * Kepler Multi-body & Transit Timing Variations Working Groups
  * Origins of Solar Systems program, award NNX14AI76G
  * Exoplanets Research Program, award NNX15AE21G
* The Pennsylvania State University
  * Dept. of Astronomy & Astrophysics
  * Center for Exoplanets & Habitable Worlds
  * Eberly College of Science
  * Institute for CyberScience
  * Center for Astrostatistics
  * Penn State Astrobiology Research Center
* Florida Institute of Technology
* University of Florida
* Statistical and Applied Mathematical Sciences Institute
