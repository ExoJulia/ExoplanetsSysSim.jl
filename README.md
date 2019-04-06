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
* Clone the [SysSimData repository](https://github.com/ExoJulia/SysSimData). I suggest somewhere outside of your home directory, .jula  or JULIA_DEPOT_PATH. It "should" automatically download large files via git-lfs.  If not, then you can cd into the directory and run `git lfs fetch` to force it to update.
```sh
git clone git@github.com:ExoJulia/SysSimData.git 
cd SysSimData
git lfs fetch # if the data hasn't been populated
```
* Make sure that your JULIA_DEPOT_PATH (~/.julia by default) does not include an old version of CORBITS or ExopalnetsSysSim.  If this is your first time using julia v1.0, then you probably don't need to do anything.  Otherwise, I see two ways to do this:
   - One way to avoid conflicts is to move or delete the JULIA_DEPOT_PATH.  But if there's _any chance_ that you might have things in your current CORBITS or ExoplanetsSysSim repots that you want to keep, then make a backup copy of thsoe repos before deleting them.  Simillarly, if there are any other packages you've been developing, make a backup copy before deleting your JULIA_DEPOT_PATH.  Once you've fully cleared out the old repos, then 'Pkg.rm("CORBITS"); Pkg.rm("ExoplanetsSysSim"); Pkg.gc()' and rm -rf both from the dev subdirectory of your JULIA_DEPOT_PATH (~/.julia by default).  Warning:  Sometimes Julia manages to keep these around despite my efforts to delete them, so I've found it's easier to rename my .julia directory and then copy any other repos in development mode back to my new .julia directory.
   - Another way to avoid conflicts with old versions is to sepcify a new JULIA_DEPOT_PATH.  However, if you go this route, then you'll need to make sure that this environment variable is set to the desired depot in each of your future shell sessions. 
```sh
export JULIA_DEPOT_PATH=~/.julia_clean
```
One ICS-ACI, it's useful to set your JULIA_DEPOT_PATH to be in your work directory, as that is higher performance and has more space than your home directory.  E.g., 
```sh
export JULIA_DEPOT_PATH=~/work/.julia
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
* Create a symlink so data points to the SysSimData repo. 
   - Change into the directory where you're developing ExoplanetSysSim (${JULIA_DEPOT_PATH}/dev/ExoplanetsSysSim).  
   - Create a symlink named data 
```sh
cd .julia/dev/ExoplanetsSysSim
#cd ${JULIA_DEPOT_PATH}/dev/ExoplanetsSysSim  # if you set JULIA_DEPOT_PATH
ln -s PATH_TO_SYSSIMDATA data
```
   - Alternatively, you can override the default file paths to point to whereever you placed the binary input files.  Although this probably require more work. 

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
