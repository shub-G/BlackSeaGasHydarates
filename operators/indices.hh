class Indices {

public:

	const static int numOfPVs 	= 7;
	const static int PVId_Pw	= 0;
	const static int PVId_Sg	= 1;
	const static int PVId_XCH4	= 2;
	const static int PVId_YH2O	= 3;
	const static int PVId_Xc	= 4;
	const static int PVId_Sh	= 5;
	const static int PVId_T 	= 6;

	const static int numOfComps  = 4;
	const static int compId_XCH4 = 0;
	const static int compId_XH2O = 1;
	const static int compId_YCH4 = 2;
	const static int compId_YH2O = 3;

	const static int numOfSVs 	= 17;
	const static int SVId_Pg	= 0;  // gas phase pressure
	const static int SVId_Pw	= 1;  // water phase pressure
	const static int SVId_Pc	= 2;  // capillary pressure
	const static int SVId_Sg	= 3;  // gas saturation
	const static int SVId_Sw	= 4;  // water saturation
	const static int SVId_Sh	= 5;  // hydrate saturation
	const static int SVId_XCH4	= 6;  // CH4 mole fraction in water
	const static int SVId_XH2O	= 7;  // H2O mole fraction in water
	const static int SVId_YCH4	= 8;  // CH4 mole fraction in gas
	const static int SVId_YH2O	= 9;  // H2O mole fraction in gas
	const static int SVId_Xc	= 10; // salt mole fraction in water
	const static int SVId_T		= 11; // temperature
 	const static int SVId_K 	= 12; // absolute permeability
	const static int SVId_z		= 13; // methane gas compressibility factor
	const static int SVId_Peq 	= 14; // hydrate equilibrium pressure
	const static int SVId_por	= 15; // total porosity

	const static int numOfBCs = 4;
	const static int BCId_water = 0;
	const static int BCId_gas	= 1;
	const static int BCId_salt 	= 2;
	const static int BCId_heat	= 3;
	const static int neumann	= 0;
	const static int dirichlet	= 1;

#ifdef PLOT_VELOCITIES
	const static int numOfPhases = 2;
	const static int phaseId_g = 0;
	const static int phaseId_w = 1;
#endif

};
