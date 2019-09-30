# ExoplanetsSysSim
[![DOI](https://zenodo.org/badge/179777476.svg)](https://zenodo.org/badge/latestdoi/179777476)

Welcome to the ExoplanetsSysSim package for generating planetary systems and simulating observations of those systems with a transit survey.  Currently, SysSim focuses on NASA's Kepler mission, but we've aimed to develop a generic framework that can be applied to other surveys (e.g., K2, TESS, PLATO, LSST, etc.).

## How to install ExoplanetsSysSim:
* Make sure you have julia (v0.7 or greater) installed.  It's been tested on Julia v1.1.0.  You can download Julia [here](https://julialang.org/downloads/).

* Make sure you have a recent git and [git-lfs](https://git-lfs.github.com/) installed.
If you're using ICS-ACI, then do this by running the following for each shell (or putting it in your .bashrc)
```sh
export PATH=/gpfs/group/dml129/default/sw/git-lfs:$PATH
module load git
```
* If you want to use ssh keys instead of https authentication (to minimize typing your github password), then:
  * Setup a local ssh key using ssh-keygen
  * Tell Github about your ssh key:  Person Icon (upper right), Settings, SSH & GPG keys, New SSH Key.  Enter a name in the title box and paste the contents of `cat ~/.ssh/id_rsa.pub` into the "Key" box. Add SSH Key.

* Create a clone of the [SysSimData repository](https://github.com/ExoJulia/SysSimData). 
   - If you might want to add/edit files in the SysSimData repository, then please fork your own repository on github and clone that instead of the repo in ExoJulia.  Then you can create pull requests when you're ready to add/update files in the main repository.  
   - If you plan to use existing SysSimData files, then you can just create a new copy, use `git clone`.  I suggest somewhere outside of your home directory, .julia  or JULIA_DEPOT_PATH.  
Once you've got a clone of a SysSimData repository, initialize and update the submodules.  Git "should" automatically download large files via git-lfs.  If not, then you can cd into the directory and run `git lfs fetch` to force it to update.  For example, 
```sh
git clone git@github.com:ExoJulia/SysSimData.git 
cd SysSimData
git submodule init
git submodule update
git lfs fetch # if the binary data files didn't download automatically
```
   - If you're using ICS-ACI, then you could simply use the repo in /storage/home/ebf11/group/ebf11/kepler/SysSimData that should already be set up

* Make sure that your JULIA_DEPOT_PATH (~/.julia by default) does not include an old version of CORBITS or ExopalnetsSysSim.  If this is your first time using julia v1.0, then you probably don't need to do anything.  Otherwise, I see two ways to do this:
   - One way to avoid conflicts is to move or delete the JULIA_DEPOT_PATH.  But if there's _any chance_ that you might have things in your current CORBITS or ExoplanetsSysSim repots that you want to keep, then move rather than delete (or make a backup copy of those repos before deleting them).  Simillarly, if there are any other packages you've been developing, make sure you have a backup copy before deleting your JULIA_DEPOT_PATH.            Once you've fully cleared out the old repos, then 'Pkg.rm("CORBITS"); Pkg.rm("ExoplanetsSysSim"); Pkg.gc()' and 'rm -rf CORBITS ExoplanetsSysSim' both from the dev subdirectory of your JULIA_DEPOT_PATH (~/.julia by default).  Warning:  Sometimes Julia manages to keep these around despite my efforts to delete them, so I've found it's easier to rename my .julia directory and then copy any other repos in development mode back to my new .julia directory.
   - Another way to avoid conflicts with old versions is to sepcify a new JULIA_DEPOT_PATH.  However, if you go this route, then you'll need to make sure that this environment variable is set to the desired depot in each of your future shell sessions. 
```sh
export JULIA_DEPOT_PATH=~/.julia_clean
```
On ICS-ACI, it's useful to set your JULIA_DEPOT_PATH to be in your work directory, as that is higher performance and has more space than your home directory.  I've put this in my .bashrc, so I don't forget and get confused about what's being modified.  E.g., 
```sh
export JULIA_DEPOT_PATH=~/work/.julia
```

* Run julia and install the ExoplanetsSysSim repo as a Julia package.  
  - If you will only be using it as is, then you can simply add the registered repo under the ExoJulia organization.
```julia
using Pkg
Pkg.add("ExoplanetsSysSim")
```
  - However, if you may be modifying source code in the ExoplanetsSysSim directory itself, then please fork your own version on github and develop that version instead.  For example,
```julia
Pkg.develop(PackageSpec(url="git@github.com:ExoJulia/ExoplanetsSysSim.jl.git"))
```
(but replacing ExoJulia with the github username associated with your fork).  If you've set ExoplanetsSysSim to be under development, Julia will not automatically update it.  You'll have to do a `git pull` from dev/ExoplanetsSysSim to merge in new updates.

  - Some MacOS users find that CORBITS does not build successfully.  This does not prevent MacOS users from using SysSim in "single-observer mode" (which is the mode used for existing publications).

* Create a symlink so 'data' in the ExoplanetsSysSim directory points to the SysSimData repo. 
   - Change into the directory where you've added or developing ExoplanetSysSim (likely ${JULIA_DEPOT_PATH}/dev/ExoplanetsSysSim).  
   - Create a symlink named data 
```sh
cd .julia/dev/ExoplanetsSysSim
#cd ${JULIA_DEPOT_PATH}/dev/ExoplanetsSysSim  # alternative if you set JULIA_DEPOT_PATH
ln -s PATH_TO_SYSSIMDATA data
```
   - Alternatively, you can override the default file paths to point to wherever you placed the binary input files.  Although this probably require more work. 

* Optionally, run some tests, e.g.
```julia
using ExoplanetsSysSim
include(joinpath(dirname(pathof(ExoplanetsSysSim)),"..","test","runtests.jl"))
```
## How to use SysSim for your own Projects
- Install ExoplanetsSysSim (see above)
- Create your own repository containing code that will call ExoplanetsSysSim
- Make it a Julia project by adding dependencies, including ExoplanetsSysSim.
- Make your project depend on either the registered version of ExoplanetsSysSim or the version in your development directory.  Since you've already installed ExoplanetSysSim, then Julia should find and reuse the code in the dev directory rather than reinstalling it. 
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.add("ExoplanetsSysSim")  # For the registered version of ExoplanetsSysSim
# Pkg.develop("ExoplanetsSysSim") # To use your development branch of ExoplanetsSysSim.
```   
- Have your project code load ExoplanetsSysSim and use it
```julia
using ExoplanetsSysSim
...
```
   - At the moment, you can test using 'examples/generate_catalogs.jl' from Matthias's project at https://github.com/ExoJulia/SysSimExClusters
   - By default, the master branch includes recent updates.  There is a chance that we occasionally break something that's not part of our test set.  Therefore, we've created a [stable branch](https://github.com/ExoJulia/ExoplanetsSysSim.jl/tree/stable) which users may wish to use for science results to be published.  If you find something broken in the stable branch, then please check the [open issues](https://github.com/ExoJulia/ExoplanetsSysSim.jl/issues).  If we're not yet aware of your problem, then notify the SysSim team via a new GitHub issue.
   
* Write your papers and share your code as a GitHub repo
   - If you want to share your Manifest.toml file, then make a copy of the Manifest.toml when you're not in develop mode.  Otherwise, users on other systems will get errors, since they can't access the same path with your development version.
   - If you'd like your code to appear as part of the [ExoJulia organization](https://github.com/ExoJulia/), then please let [Eric](https://github.com/eford) know.

* Cite relevant code and associated publications
  - [![DOI](https://zenodo.org/badge/179777476.svg)](https://zenodo.org/badge/latestdoi/179777476) Citation for core SysSim code base.
  - [Hsu et al. (2018) AJ 155, 205.](https://arxiv.org/ct?url=https%3A%2F%2Fdx.doi.org%2F10.3847%2F1538-3881%2Faab9a8&v=19ae32f8) (first published paper, describes basic SysSim functionality pre-1.0 version, please cite until Hsu et al. 2019 is accepted)
  - [Hsu et al. (2019) submitted to AJ. arXiv:1902.01417](https://arxiv.org/abs/1902.01417) (most recent public paper, describing improvements to model for Kepler pipeline in SysSim v1.0, please cite if using SysSim v1.*)
  - [He, Ford, & Ragozzine (2019) submitted to MNRAS. arXiv:1907.07773](https://arxiv.org/abs/1907.07773) (describes model for generating planetary systems, uses SysSim v1.0, please cite if using clustered model)
  - [Brakensiek & Ragozzine (2016) ApJ 821, 47.](https://doi.org/10.3847/0004-637X/821/1/47) and (citation for CORBITS, please cite if you make use of averaging over viewing geometries)
* Let the SysSim team know about your publication (or other use of SysSim, e.g., proposals) via pull request

## The SysSim Team:
### Key Developers:
  * Eric Ford:  Conceptual framework, Development of core codebase
  * Matthias He:  Development and application of clustered multi-planet model
  * Danley Hsu:  Validation of Kepler model, distance functions and application to planet occurence rates
  * Darin Ragozzine:  Conceptual framework, Incorporating DR25 data products
### Other Contributors/Consultants:
  * Robert Morehead:  Preliminary model development, exploratory applications of ABC and comparing distance functions.
  * Keir Ashby:  Testing incorporation of DR 25 data products
  * Jessi Cisewski:  Advice on statistical methodlogy
  * Chad Schafer:  Advice on statistical methodlogy
  * Tom Loredo:  Advice on statistical methodlogy
  * Robert Wolpert:  Advice on statistical methodlogy

### Acknowledgements:
* NASA
  * [Kepler Mission](https://www.nasa.gov/mission_pages/kepler/main/index.html)
  * [Kepler Science Team](https://www.nasa.gov/mission_pages/kepler/team/teamroster)
  * Kepler Multi-body & Transit Timing Variations Working Groups
  * Origins of Solar Systems program, award NNX14AI76G
  * Exoplanets Research Program, award NNX15AE21G
* [The Pennsylvania State University](https://www.psu.edu/)
  * [Dept. of Astronomy & Astrophysics](http://astro.psu.edu/)
  * [Center for Exoplanets & Habitable Worlds](https://exoplanets.psu.edu/)
  * [Eberly College of Science](http://science.psu.edu/)
  * [Institute for CyberScience](https://ics.psu.edu/)
  * [Center for Astrostatistics](https://astrostatistics.psu.edu/)
  * [Penn State Astrobiology Research Center](http://psarc.weebly.com/)
* [Brigham Young University](https://www.physics.byu.edu/)
* [University of Florida](https://www.ufl.edu/)
* [Florida Institute of Technology](https://www.fit.edu/)
* [Statistical and Applied Mathematical Sciences Institute](https://www.samsi.info/)
