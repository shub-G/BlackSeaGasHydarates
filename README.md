# BlackSeaGasHydarates
This repository contains the source code used in the numerical simulations presented in the manuscript "Simulating the development of multiple BSRs in the dynamic gas hydrate system of the Danube paleo-delta, Black Sea, over the past 300 ka.", authored by S. Gupta, C. Deusner, E.B.-Galerne, and M. Haeckel.

## Instructions for installation
* System requirements: 
  * Ububtu LTS 18.04 or higher
  * packages: git, auto-conf, cmake (>=3.13), clang, g++, gcc, gfortran, superlu (dev), lapack, libblas, libboost, metis (dev), parmetis (dev), gmsh, paraview, openmpi 
* Execute following commands in terminal (in given order):
  * mkdir dune_2_8
  * cd mkdir dune_2_8
  * mkdir source
  * cd mkdir source
  * cp /downloads/installer.sh .
  * chmod 755 installer.sh
  * ./installer.sh dune
  * cd dune
  * ./buildmobules.sh
  * cd ../..
  * chmod 755 newproject.sh
  * ./newproject.sh BlackSeaGasHydrates
    * On prompt for "2) Which modules should this module depend on?", enter: dune-common dune-geometry dune-uggrid dune-grid dune-localfunctions dune-istl dune-typetree dune-functions dune-alugrid dune-pdelab
    * On prompt for "3) Project/Module version?", enter 1
    * On promt for "4) Maintainer's email address?", enter your email address.
  * chmod 755 buildproject.sh
  * ./buildproject.sh BlackSeaGasHydrates
  * cd BlackSeaGasHydrates/src/
  * rm -rf blacksea_multiple_bsrs.cc
  * rm -rf CMakeLists.txt
  * cp \_all_source_files_in_repo\_ .
  * cd ../..
  * chmod 755 compile.sh
  * ./compile.sh BlackSeaGasHydrates

## To run the simulations:
* Execute following commands in terminal (in given order):
  * cd /dune_2_8/BlackSeaGasHydrates/release-build/src
  * ./blacksea_multiple_bsrs \_your-user-name\_ \_input-file\_  
    * input files are located in the folder: /dune_2_8/BlackSeaGasHydrates/src/blacksea_multiple_bsrs/inputs/
    * input files are the files with the extension ".ini". In the execution call, drop the .ini extension.
    * hint on \_your_user_name\_: The main executable looks for the following path: home/\_your_user_name\_/dune_2_8/BlackSeaGasHydrates/

## Files included in this repo:
* installation files:
  * installation/installer.sh
  * installation/newproject.sh
  * installation/buildproject.sh
  * installation/compile.sh
* source files 
  * include_dune.hh
  * blacksea_multiple_bsrs.cc
  * blacksea_multiple_bsrs/include_problem.hh
    * This includes all the necessary src files
  * blacksea_multiple_bsrs/inputs/scenarioX.ini 
    * Input files are set up for scenarios with each combination of the parameters sampled in the study presented in this manuscript.
    * These parameters are: 1) initial GH volume V0, 2) kinetic rate of hydrate dissociation kd, and 3) kinetic rate of hydrate formation kf.
* outputs/blacksea_multiple_bsrs/test0
  *  Simulation outputs for the reference scenario (scenario2) are archived in this repository.
  *  Use PARAVIEW to visualize and plot this scenario.
* Animations of the reference scenario
  *  outputs/blacksea_multiple_bsrs/BlackSea_multipleBSRs_scenario2_small_fast.avi
  *  outputs/blacksea_multiple_bsrs/BlackSea_multipleBSRs_scenario2_small_fast.mp4
