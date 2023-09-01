/*
 * blacksea_multiple_bsrs.cc
 *
 *  Created on: Apr 12, 2022
 *      Author: sgupta
 */

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>

/*********************************************************************/
/* SETTINGS
 *
 * options:
 * (grig backend -> USE_YSP, USE_UG)
 * (mpi -> PARALLEL, SEQUENTIAL)
 * (properties -> FRAUENHOFER_MODEL, PENG_ROBINSON_EOS)
 * (problem IDs: ...)
 *
 * default options: USE_YASP, SEQUENTIAL, BASE_EOS
 */
#define PARALLEL
#define USE_UG
#define FRAUENHOFER_MODEL
#define COMPACTION
#define PLOT_VELOCITIES

/* PROBLEM DESCRIPTION
 * 		2D problem setting with the Channel Levee geometry based on the geological setting of the Black Sea.
 * 		Sedimentation induced rise in sea floor temperature resulting in hydrate dissociation and gas migration through GHSZ.
 * for sequential run:
 * ./blacksea_multiple_bsrs <user_name> <inputs_file>
 */

/*********************************************************************/
#include"include_dune.hh"
#include"blacksea_multiple_bsrs/include_problem.hh"
/*********************************************************************/

int main(int argc, char** argv)
{
	  try{
		// Maybe initialize MPI
		Dune::MPIHelper& helper = Dune::MPIHelper::instance(argc, argv);
		if(helper.rank()==0){
			std::cout << "Hello World! This is problem blacksea_multiple_bsrs of HydrateProject." << std::endl;
		}
		if(Dune::MPIHelper::isFake){
			std::cout<< "This is a sequential program." << std::endl;
		}
		else {
			if(helper.rank()==0){
				std::cout<<"I am rank "<<helper.rank()<<" of "<<helper.size()<<" processes!"<<std::endl;
			}
		}

		/**************************************************************************************************/
		// INPUTS
	    	if (argc!=2)
	    	{
		    	if(helper.rank()==0){
		    		std::cout << "usage: ./blacksea_multiple_bsrs <user_name> <inputs_file.ini>" << std::endl;
		    	}
			return 1;
	    	}
	        /**************************************************************************************************/
		// USER-NAME
		char USER_NAME[100];
		sscanf(argv[1],"%99s", USER_NAME);
		// DUNE MODEL PATH
		std::string MODEL_PATH = "/home/";
		MODEL_PATH += USER_NAME;
		MODEL_PATH += "/dune_2_8/HydrateProjectDemo/src/";
		// PROBLEM NAME
		std::string PROBLEM_NAME = "blacksea_multiple_bsrs/";
		// INPUT PATH NAME
		std::string INPUT_PATH  = MODEL_PATH + PROBLEM_NAME + "inputs/";
		// OUTPUT PATH NAME
		std::string OUTPUT_PATH = MODEL_PATH + "outputs/" + PROBLEM_NAME;
		// INI-FILE FOR USER-DEFINED INPUTS
		char INI_FILE[100];
		sscanf(argv[2],"%99s", INI_FILE);
		std::string input_file = INPUT_PATH;
		input_file += INI_FILE;
		input_file += ".ini";
		if(helper.rank()==0){
		std::cout<< "input file: " << input_file << std::endl ;
		}
		/**************************************************************************************************/
		// PARAMETER TREE
		Dune::ParameterTree ptree;
		Dune::ParameterTreeParser ptreeparser;
		ptreeparser.readINITree(input_file,ptree);
		ptreeparser.readOptions(argc,argv,ptree);

		/**************************************************************************************************/
		// MESH
		typedef std::vector<int> GmshIndexMap;
		GmshIndexMap boundary_index_map;
		GmshIndexMap element_index_map;

		// MESH
	    	MeshParameters<Dune::ParameterTree> mesh(ptree);
	    	const int dim = mesh.dimension;

#ifdef USE_YASP
		/*************
		 *  YASP
		 *************/
	    	Dune::FieldVector<double,dim> L;
        	L[0] = mesh.X_length; /*ndim*/
        	L[1] = mesh.Z_length; /*ndim*/
        	std::array<int,dim> N;
        	N[0] = mesh.X_cells;
       	N[1] = mesh.Z_cells;

		std::bitset<dim> periodic(false);
		int overlap=1;
        	typedef Dune::YaspGrid<dim> Grid;
        	std::shared_ptr<Grid> grid = std::shared_ptr<Grid>(new Grid(L,N,periodic,overlap,helper.getCommunicator()));
        	typedef Grid::LeafGridView GV;
       	GV gv=grid->leafGridView();
        	grid->loadBalance();

#elif defined(USE_UG)
		/*************
		 *  UG
		 *************/
		typedef Dune::UGGrid<dim> GridType;
		GridType grid_type;
		const std::string grid_name = ptree.get("grid.ug.name", (std::string)"grid") ;
		std::string grid_file = MODEL_PATH + PROBLEM_NAME + "grids/";
		grid_file += grid_name;
		grid_file += ".msh";
		Dune::GmshReader<GridType> gmshreader;
		std::shared_ptr<GridType> grid(gmshreader.read(grid_file,boundary_index_map, element_index_map,true,false));

		typedef GridType::LeafGridView GV;
		GV gv = grid->leafGridView();
        	grid->loadBalance();

#else
        	std::cout<< "Incorrect grid manager. Please choose YASP or UG." << std::endl;
        	exit(0);
#endif

		/**************************************************************************************************/
		// DRIVER
		driver( gv, ptree,
			boundary_index_map,
			element_index_map,
        		OUTPUT_PATH,
			helper);

		/**************************************************************************************************/

	  }
	  catch (Dune::Exception &e){
	    std::cerr << "Dune reported error: " << e << std::endl;
	  }
	  catch (...){
	    std::cerr << "Unknown exception thrown!" << std::endl;
	  }
}
