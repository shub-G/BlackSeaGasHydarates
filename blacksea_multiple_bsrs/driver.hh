template< typename GV, typename PTree, typename GmshIndexMap >
void driver( 	const GV& gv, 				// GridView
		const PTree& ptree, 			// Input Parameters
		GmshIndexMap boundary_index_map,	// Boundary masks
		GmshIndexMap element_index_map,	// Element masks
		std::string output_path,
		Dune::MPIHelper& helper){

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
		typedef typename GV::Grid::ctype Coord;
		const int dim = GV::dimension;

		double time = 0.0;
		double dt = 0.0;

	//	MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
		typedef Properties<GV,PTree,GmshIndexMap> Properties;
		Properties property(gv,ptree,element_index_map,&time,&dt);
		
		std::string pathName = output_path;
		auto pathExt = ptree.get("output.path_name",(std::string)"test0");
		pathName += pathExt;
		pathName += "/";
		auto fileName = ptree.get("output.file_name",(std::string)"test");

		/*Non-dimensionalize time prams*/
		/*in input file dt, t_END, dt_min, dt_max, t_OP are specified in kilo-annum*/
		dt  = ptree.get("time.dt_initial",(double)0.001);
		dt *= (1000.*365.*24.*60.*60.); /*convert to seconds*/
		dt *= 1./property.characteristicValue.t_c; /*ndim*/
		double t_END  = ptree.get("time.time_end",(double)300.0);
		t_END *= (1000.*365.*24.*60.*60.); /*convert to seconds*/
		t_END *= 1./property.characteristicValue.t_c; /*ndim*/
		// output time interval
		double t_OP   = ptree.get("output.time_interval",(double)1000);
		t_OP *= (1000.*365.*24.*60.*60.); /*convert to seconds*/
		t_OP *= 1./property.characteristicValue.t_c; /*ndim*/
		//adaptive time control
		bool adaptive_time_control = ptree.get("adaptive_time_control.flag",(bool)true);
		double dt_min = ptree.get("adaptive_time_control.dt_min",(double)1.e-6);
		dt_min *= (1000.*365.*24.*60.*60.); /*convert to seconds*/
		dt_min *= 1./property.characteristicValue.t_c; /*ndim*/
		double dt_max = ptree.get("adaptive_time_control.dt_max",(double)10.);
		dt_max *= (1000.*365.*24.*60.*60.); /*convert to seconds*/
		dt_max *= 1./property.characteristicValue.t_c; /*ndim*/
		int maxAllowableIterations = ptree.get("adaptive_time_control.max_newton_steps",(int)6);
		int minAllowableIterations = ptree.get("adaptive_time_control.min_newton_steps",(int)4);

		double dtstart = dt;
		double time_op = time;
		double clock_time_elapsed = 0.;

		/************************************************************************************************/
		// MAIN
		/************************************************************************************************/
		//	COMPOSITE GFS FOR PRIMARY VARIABLES
		typedef typename GV::Grid::ctype Coord;
#ifdef PARALLEL
		typedef Dune::PDELab::P0ParallelConstraints CON0;
#else
		typedef Dune::PDELab::NoConstraints CON0;
#endif
		using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;					// default block size: 1
		typedef Dune::PDELab::QkDGLocalFiniteElementMap<Coord,double,0,dim,Dune::PDELab::QkDGBasisPolynomial::lagrange> FEM0;
		FEM0 fem0;
		typedef Dune::PDELab::GridFunctionSpace<GV, FEM0, CON0, VBE0> GFS0;
		GFS0 gfs0(gv, fem0);
		// gfs for composite system: Pw , Sg , XCH4 , YH2O , Xc , Sh , T
#ifdef PARALLEL
		using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  Indices::numOfPVs,
													  VBE,
													  Dune::PDELab::EntityBlockedOrderingTag > GFS;
		GFS gfs(gfs0);
#else
		using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  Indices::numOfPVs,
													  VBE,
													  Dune::PDELab::LexicographicOrderingTag > GFS;
		GFS gfs(gfs0);
#endif
		typedef typename GFS::template ConstraintsContainer<double>::Type CC;
		CC cc;
		cc.clear();

		//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
		using PathPw = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Pw>>;
		using SUBGFS_Pw = Dune::PDELab::GridFunctionSubSpace<GFS,PathPw>;
		SUBGFS_Pw	subgfs_Pw(gfs);
		using PathSg = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sg>>;
		using SUBGFS_Sg = Dune::PDELab::GridFunctionSubSpace<GFS,PathSg>;
		SUBGFS_Sg	subgfs_Sg(gfs);
		using PathXCH4 = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_XCH4>>;
		using SUBGFS_XCH4 = Dune::PDELab::GridFunctionSubSpace<GFS,PathXCH4>;
		SUBGFS_XCH4	subgfs_XCH4(gfs);
		using PathYH2O = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_YH2O>>;
		using SUBGFS_YH2O = Dune::PDELab::GridFunctionSubSpace<GFS,PathYH2O>;
		SUBGFS_YH2O	subgfs_YH2O(gfs);
		using PathXc = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Xc>>;
		using SUBGFS_Xc = Dune::PDELab::GridFunctionSubSpace<GFS,PathXc>;
		SUBGFS_Xc	subgfs_Xc(gfs);
		using PathSh = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_Sh>>;
		using SUBGFS_Sh = Dune::PDELab::GridFunctionSubSpace<GFS,PathSh>;
		SUBGFS_Sh	subgfs_Sh(gfs);
		using PathT = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_T>>;
		using SUBGFS_T = Dune::PDELab::GridFunctionSubSpace<GFS,PathT>;
		SUBGFS_T	subgfs_T(gfs);

		//	MAKE VECTOR CONTAINER FOR THE SOLUTION
		using U = Dune::PDELab::Backend::Vector<GFS,double>;
		U u_old(gfs);
		U u_new(gfs);

		//  EVALUATION FUNCTIONS FOR PRIMARY VARIABLES
		Dune::PDELab::Evaluation<SUBGFS_Pw	,U> evaluation_Pw(	 subgfs_Pw	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_Sg	,U> evaluation_Sg(	 subgfs_Sg	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_XCH4,U> evaluation_XCH4( subgfs_XCH4,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_YH2O,U> evaluation_YH2O( subgfs_YH2O,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_Xc	,U> evaluation_Xc(	 subgfs_Xc	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_Sh	,U> evaluation_Sh(	 subgfs_Sh	,&u_new );
		Dune::PDELab::Evaluation<SUBGFS_T	,U> evaluation_T( 	 subgfs_T 	,&u_new );

		//	INITIAL CONDITIONS
		//	Make function for initial values
		typedef Pw_Initial<GV,Properties,GmshIndexMap> ICV_Pw;
		ICV_Pw Pw_initial(gv,property,element_index_map);
		typedef Sg_Initial<GV,Properties,GmshIndexMap> ICV_Sg;
		ICV_Sg Sg_initial(gv,property,element_index_map);
		typedef XCH4_Initial<GV,Properties,GmshIndexMap> ICV_XCH4;
		ICV_XCH4 XCH4_initial(gv,property,element_index_map);
		typedef YH2O_Initial<GV,Properties,GmshIndexMap> ICV_YH2O;
		ICV_YH2O YH2O_initial(gv,property,element_index_map);
		typedef Xc_Initial<GV,Properties,GmshIndexMap> ICV_Xc;
		ICV_Xc Xc_initial(gv,property,element_index_map);
		typedef Sh_Initial<GV,Properties,GmshIndexMap> ICV_Sh;
		ICV_Sh Sh_initial(gv,property,element_index_map);
		typedef T_Initial<GV,Properties,GmshIndexMap> ICV_T ;
		ICV_T T_initial( gv,property,element_index_map);
		typedef Dune::PDELab::CompositeGridFunction< ICV_Pw,
													 ICV_Sg,
													 ICV_XCH4,
													 ICV_YH2O,
													 ICV_Xc,
													 ICV_Sh,
													 ICV_T  > InitialValues;
		InitialValues icv( Pw_initial, Sg_initial, XCH4_initial, YH2O_initial, Xc_initial, Sh_initial, T_initial );

		// 	Initialize the solution at t=0 (uold) with the given initial values
		Dune::PDELab::interpolate( icv, gfs, u_old );
		u_new = u_old;

		//	BOUNDARY CONDITIONS
		typedef ProblemBoundaryConditions<GV,Properties,GmshIndexMap> BoundaryConditions ;
		BoundaryConditions bc( gv,property,boundary_index_map ) ;

		//	MAKE INSTATIONARY GRID OPERATOR SPACE

		typedef LocalOperator< GV, Properties, BoundaryConditions > LOP;	// spatial part
		LOP lop( gv, property, bc );
		typedef TimeOperator< GV, Properties > TLOP; // temporal part
		TLOP tlop( gv, property );

		typedef Dune::PDELab::ISTL::BCRSMatrixBackend<> MBE;
		MBE mbe(35);

		typedef Dune::PDELab::GridOperator< GFS, GFS, LOP, MBE, double, double, double, CC, CC> GOLOP;
		GOLOP golop(gfs, cc, gfs, cc, lop, mbe);

		// How well did we estimate the number of entries per matrix row?
		// => print Jacobian pattern statistics
		typename GOLOP::Traits::Jacobian jac(golop);
		if(helper.rank()==0){
			std::cout << jac.patternStatistics() << std::endl;
		}

		typedef Dune::PDELab::GridOperator< GFS, GFS, TLOP, MBE, double, double, double, CC, CC > GOTLOP;
		GOTLOP gotlop(gfs, cc, gfs, cc, tlop, mbe);

		typedef Dune::PDELab::OneStepGridOperator< GOLOP, GOTLOP > IGO;
		IGO igo( golop, gotlop );

		// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL
		typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO> LS; //works
		LS ls(gfs,100,1,false,true);
#else
		typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LS; //works
		LS ls(0);
#endif

		//	SELECT SOLVER FOR NON-LINEAR PROBLEM
		typedef Dune::PDELab::Newton< IGO, LS, U > PDESOLVER;
		PDESOLVER pdesolver( igo, ls );
		// 	select control parameters for non-linear PDE-solver
		pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::noLineSearch);
		pdesolver.setReassembleThreshold(0.0);
		pdesolver.setVerbosityLevel(2);
		pdesolver.setReduction(ptree.get("newton.reduction",(double)1.e-6));
		pdesolver.setMinLinearReduction(1e-9);
		pdesolver.setMaxIterations(ptree.get("newton.max_iterations",(int)15));
		pdesolver.setForceIteration(true);
		pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error",(double)1.e-6));


		//	SELECT TIME-STEPPER
		Dune::PDELab::ImplicitEulerParameter<double> method;

		Dune::PDELab::OneStepMethod< double, IGO, PDESOLVER, U, U > osm( method, igo, pdesolver );
		osm.setVerbosityLevel(2);

		/************************************************************************************************/
		//  POST-PROCESS - PV/SV: Pg, Pw, Pc, Sg, Sw, Sh, XCH4, XH2O, YCH4, YH2O, Xc, T, K, zCH4, Peq, por
		/************************************************************************************************/
		typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none > VBE_PP;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  Indices::numOfSVs,
													  VBE_PP,
													  Dune::PDELab::LexicographicOrderingTag > GFS_PP;
		GFS_PP gfs_pp(gfs0);
		typedef typename GFS_PP::template ConstraintsContainer<double>::Type CC_PP;
		CC_PP cc_pp;
		cc_pp.clear();

		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pg	>> > SUBPP_Pg;
		SUBPP_Pg 	subpp_Pg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pw	>> > SUBPP_Pw;
		SUBPP_Pw 	subpp_Pw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Pc	>> > SUBPP_Pc;
		SUBPP_Pc 	subpp_Pc(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sg	>> > SUBPP_Sg;
		SUBPP_Sg 	subpp_Sg(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sw	>> > SUBPP_Sw;
		SUBPP_Sw 	subpp_Sw(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Sh	>> > SUBPP_Sh;
		SUBPP_Sh 	subpp_Sh(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XCH4	>> > SUBPP_XCH4;
		SUBPP_XCH4 	subpp_XCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_XH2O >> > SUBPP_XH2O;
		SUBPP_XH2O 	subpp_XH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YCH4 >> > SUBPP_YCH4;
		SUBPP_YCH4 	subpp_YCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_YH2O >> > SUBPP_YH2O;
		SUBPP_YH2O 	subpp_YH2O(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Xc   >> > SUBPP_Xc;
		SUBPP_Xc 	subpp_Xc(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_T	>> > SUBPP_T;
		SUBPP_T 	subpp_T(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_K	>> > SUBPP_K;
		SUBPP_K 	subpp_K(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_z	>> > SUBPP_zCH4;
		SUBPP_zCH4 	subpp_zCH4(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_Peq	>> > SUBPP_Peq;
		SUBPP_Peq 	subpp_Peq(gfs_pp);
		typedef typename Dune::PDELab::GridFunctionSubSpace< GFS_PP, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::SVId_por	>> > SUBPP_por;
		SUBPP_por 	subpp_por(gfs_pp);

		using U_PP = Dune::PDELab::Backend::Vector<GFS_PP,double>;
		U_PP u_pp(gfs_pp, 0.0);

		PostProcess< GV, Properties,
					 Dune::PDELab::Evaluation<SUBGFS_Pw  ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_Sg	 ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_XCH4,U>,
					 Dune::PDELab::Evaluation<SUBGFS_YH2O,U>,
					 Dune::PDELab::Evaluation<SUBGFS_Xc  ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_Sh  ,U>,
					 Dune::PDELab::Evaluation<SUBGFS_T   ,U>,
					 GFS_PP, U_PP > postprocess( gv, property,
												 &evaluation_Pw,
												 &evaluation_Sg,
												 &evaluation_XCH4,
												 &evaluation_YH2O,
												 &evaluation_Xc,
												 &evaluation_Sh,
												 &evaluation_T,
												 gfs_pp, &u_pp);
		postprocess.evaluate();

#ifdef PLOT_VELOCITIES
		/************************************************************************************************/
		//  POST-PROCESS - PHASE VELOCITIES
		//	step 1: project phase pressures onto Q1 (C1-continuous) space (LOP_PROJ, SLP_PROJ)
		//	step 2: evaluate phase velocities at cell centers (PHASE_VELOCITY)
		/************************************************************************************************/

		// STEP 1:
		// L2-PROJECTION OF P0-PRESSURES ONTO Q1-SPACE
		const int degree = 1; //polynomial degree
		typedef Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, double, degree> FEMPROJ;
		FEMPROJ femproj(gv);
		typedef Dune::PDELab::GridFunctionSpace<GV, FEMPROJ, CON0, VBE0> GFSPROJ0;
		GFSPROJ0 gfsproj0(gv, femproj);
#ifdef PARALLEL
		using VBEPROJ = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFSPROJ0,
													  Indices::numOfPhases,
													  VBEPROJ,
													  Dune::PDELab::EntityBlockedOrderingTag > GFSPROJ;
		GFSPROJ gfsproj(gfsproj0);
#else
		typedef Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none > VBEPROJ;
		typedef Dune::PDELab::PowerGridFunctionSpace< GFSPROJ0,
													  Indices::numOfPhases,
													  VBEPROJ,
													  Dune::PDELab::LexicographicOrderingTag > GFSPROJ;
		GFSPROJ gfsproj(gfsproj0);
#endif
		typedef typename GFSPROJ::template ConstraintsContainer<double>::Type CCPROJ;
		CCPROJ ccproj;
		ccproj.clear();

		using SUBPROJ_Pg = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::phaseId_g>> >;
		SUBPROJ_Pg subproj_Pg(gfsproj);
		using SUBPROJ_Pw = Dune::PDELab::GridFunctionSubSpace< GFSPROJ, Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::phaseId_w>> >;
		SUBPROJ_Pw subproj_Pw(gfsproj);

		using UPROJ = Dune::PDELab::Backend::Vector<GFSPROJ,double>;
		UPROJ u_proj(gfsproj); // u_proj = [Pg_q1, Pw_q1]^T

		Dune::PDELab::Evaluation<SUBPROJ_Pg	,UPROJ> evaluation_projPg(	 subproj_Pg	,&u_proj );
		Dune::PDELab::Evaluation<SUBPROJ_Pw	,UPROJ> evaluation_projPw(	 subproj_Pw	,&u_proj );

		typedef LocalOperatorPROJ< GV,Properties,GFS,U> LOPPROJ;	// spatial part
		LOPPROJ lopproj(gv, property, gfs, &u_new);

		MBE mbeproj(10);
		typedef Dune::PDELab::GridOperator< GFSPROJ, GFSPROJ, LOPPROJ, MBE, double, double, double, CCPROJ, CCPROJ> GOPROJ;
		GOPROJ goproj(gfsproj, ccproj, gfsproj, ccproj, lopproj, mbeproj);

#ifdef PARALLEL
	#ifdef USE_UG
			typedef Dune::PDELab::ISTLBackend_OVLP_GMRES_ILU0<GFSPROJ,CCPROJ> LSPROJ;
			LSPROJ lsproj(gfsproj,ccproj,100,1,200);
	#elif defined(USE_YASP)
			typedef Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<GOPROJ> LSPROJ;
			LSPROJ lsproj(gfsproj,100,1,true,true);
	#endif
#else
		typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LSPROJ;
		LSPROJ lsproj;
#endif

		typedef Dune::PDELab::StationaryLinearProblemSolver<GOPROJ,LSPROJ,UPROJ> SLPPROJ;
		SLPPROJ slpproj(goproj,lsproj,u_proj,1e-10);
		slpproj.apply(); // get initial values of u_proj

		// STEP 2:
		// PHASE VELOCITY CALCULATOR
		// 2.1->
		typedef Dune::PDELab::PowerGridFunctionSpace< GFS0,
													  dim,
													  Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none >,
													  Dune::PDELab::LexicographicOrderingTag > GFS_V;
		GFS_V gfs_v(gfs0);
		typedef typename GFS_V::template ConstraintsContainer<double>::Type CC_V;
		CC_V cc_v;
		cc_v.clear();

		using U_V = Dune::PDELab::Backend::Vector<GFS_V,double>;
		U_V u_vg(gfs_v, 0.0); // gas phase velocity
		U_V u_vw(gfs_v, 0.0); // water phase velocity

		PhaseVelocity< GV, Properties,
					   Dune::PDELab::Evaluation<SUBPROJ_Pg ,UPROJ>,
					   Dune::PDELab::Evaluation<SUBGFS_Pw  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Sg  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_XCH4,U>,
					   Dune::PDELab::Evaluation<SUBGFS_YH2O,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Xc  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Sh  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_T   ,U>,
					   GFS_V, U_V > gasvelocity( gv, property,
							   	   	   	   	     &evaluation_projPg,
									  	  	  	 &evaluation_Pw,
												 &evaluation_Sg,
												 &evaluation_XCH4,
												 &evaluation_YH2O,
												 &evaluation_Xc,
												 &evaluation_Sh,
												 &evaluation_T,
												 gfs_v, Indices::phaseId_g, &u_vg,
												 &time, &dt);

		PhaseVelocity< GV, Properties,
					   Dune::PDELab::Evaluation<SUBPROJ_Pw ,UPROJ>,
					   Dune::PDELab::Evaluation<SUBGFS_Pw  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Sg  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_XCH4,U>,
					   Dune::PDELab::Evaluation<SUBGFS_YH2O,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Xc  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_Sh  ,U>,
					   Dune::PDELab::Evaluation<SUBGFS_T   ,U>,
					   GFS_V, U_V > watervelocity( gv, property,
							   	   	   	   	   	   &evaluation_projPw,
									  	  	  	   &evaluation_Pw,
												   &evaluation_Sg,
												   &evaluation_XCH4,
												   &evaluation_YH2O,
												   &evaluation_Xc,
												   &evaluation_Sh,
												   &evaluation_T,
												   gfs_v, Indices::phaseId_w, &u_vw,
												   &time, &dt);
		gasvelocity.evaluate();
		watervelocity.evaluate();

#endif

		/************************************************************************************************/
		// VTK
		/************************************************************************************************/
		// Make a grid function out of it
		//	GRAPHICS FOR INITIAL GUESS
		// secondary variables
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pg, U_PP > DGF_Pg;
		DGF_Pg dgf_Pg( subpp_Pg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pw, U_PP > DGF_Pw;
		DGF_Pw dgf_Pw( subpp_Pw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Pc, U_PP > DGF_Pc;
		DGF_Pc dgf_Pc( subpp_Pc, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sg, U_PP > DGF_Sg;
		DGF_Sg dgf_Sg( subpp_Sg, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sw, U_PP > DGF_Sw;
		DGF_Sw dgf_Sw( subpp_Sw, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Sh, U_PP > DGF_Sh;
		DGF_Sh dgf_Sh( subpp_Sh, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XCH4, U_PP > DGF_XCH4;
		DGF_XCH4 dgf_XCH4( subpp_XCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_XH2O, U_PP > DGF_XH2O;
		DGF_XH2O dgf_XH2O( subpp_XH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YCH4, U_PP > DGF_YCH4;
		DGF_YCH4 dgf_YCH4( subpp_YCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_YH2O, U_PP > DGF_YH2O;
		DGF_YH2O dgf_YH2O( subpp_YH2O, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Xc, U_PP > DGF_Xc;
		DGF_Xc dgf_Xc( subpp_Xc, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_T, U_PP > DGF_T;
		DGF_T dgf_T( subpp_T, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_K, U_PP > DGF_K;
		DGF_K dgf_K( subpp_K, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_zCH4, U_PP > DGF_zCH4;
		DGF_zCH4 dgf_zCH4( subpp_zCH4, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_Peq, U_PP > DGF_Peq;
		DGF_Peq dgf_Peq( subpp_Peq, u_pp );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPP_por, U_PP > DGF_por;
		DGF_por dgf_por( subpp_por, u_pp );
#ifdef PLOT_VELOCITIES
		typedef Dune::PDELab::DiscreteGridFunction< SUBPROJ_Pg, UPROJ > DGF_projPg;
		DGF_projPg dgf_projPg( subproj_Pg, u_proj );
		typedef Dune::PDELab::DiscreteGridFunction< SUBPROJ_Pw, UPROJ > DGF_projPw;
		DGF_projPw dgf_projPw( subproj_Pw, u_proj );
		typedef Dune::PDELab::VectorDiscreteGridFunction< GFS_V, U_V > DGF_V;
		DGF_V dgf_vg( gfs_v, u_vg );
		DGF_V dgf_vw( gfs_v, u_vw );
#endif


		int subsampling = 1;
		using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
		VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
		using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
		VTKSEQUENCEWRITER vtkSequenceWriter( std::make_shared<VTKWRITER>(vtkwriter),fileName,pathName,"");

		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pg    > >( dgf_Pg    , "Pg"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pw    > >( dgf_Pw    , "Pw"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Pc    > >( dgf_Pc    , "Pc"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sg    > >( dgf_Sg    , "Sg"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sw    > >( dgf_Sw    , "Sw"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Sh    > >( dgf_Sh    , "Sh"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_XCH4  > >( dgf_XCH4  , "XCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_XH2O  > >( dgf_XH2O  , "XH2O" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_YCH4  > >( dgf_YCH4  , "YCH4" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_YH2O  > >( dgf_YH2O  , "YH2O" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Xc    > >( dgf_Xc    , "Xc"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_T     > >( dgf_T     , "T"    ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_K     > >( dgf_K     , "K"    ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_zCH4  > >( dgf_zCH4  , "z"    ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_Peq   > >( dgf_Peq   , "Peq"  ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_por   > >( dgf_por   , "porosity"  ));
#ifdef PLOT_VELOCITIES
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_projPg > >( dgf_projPg , "Pg_proj"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_projPw > >( dgf_projPw , "Pw_proj"   ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_V > >( dgf_vg , "vg" ));
		vtkSequenceWriter.addCellData( std::make_shared< Dune::PDELab::VTKGridFunctionAdapter< DGF_V > >( dgf_vw , "vw" ));
#endif

		vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

//		//Maybe we need that, let's see
//		std::string pwrite ( const std::string & name,  const std::string & path, const std::string & extendpath,
//		                         VTK::OutputType type = VTK::ascii )


		/***********************************************/
		//	BEGIN TIME LOOP
		/***********************************************/
		int opcount = 1;
		double timecount = time;
		double dtLast = dtstart;
		int dtFlag = 0;
		bool exceptionCaught = false;
		int newton_iterations = 0;

		while( time < t_END - 1e-8/CharacteristicValues::t_c ){

			if( exceptionCaught==false ){
				dt = std::max(dt,dt_min);
			}

			if(helper.rank()==0){
			std::cout<< "_____________________________________________________" <<std::endl;
			std::cout<< " current opcount = " << opcount - 1 << std::endl;
			}

			clock_t start = clock();
			try{
				if(helper.rank()==0){
				std::cout<<"****************************" << std::endl;
				std::cout<<"  CALLING osm.apply() !"	  << std::endl;
				std::cout<<"****************************" << std::endl;
				}

				osm.apply( time, dt, u_old, u_new );

				newton_iterations = osm.getPDESolver().result().iterations;

				exceptionCaught = false;

			}catch ( Dune::Exception &e ) {
				exceptionCaught = true;
				if( dt*property.characteristicValue.t_c > 1e-8 ){

					if(helper.rank()==0){
					std::cout << "Catched Error, Dune reported error: " << e << std::endl;
					}

					u_new = u_old;

					newton_iterations = 0;

					dt *= 0.5;
					continue;
				}
				else
				{
					if(helper.rank()==0){
						std::cout << "ABORTING, due to DUNE error: " << e << std::endl;
					}
					exit(0);
				}
			}
			clock_t end = clock();
			double clock_time_this_step = (double) (end-start) / CLOCKS_PER_SEC;
			clock_time_elapsed += clock_time_this_step;

			if(helper.rank()==0){
			std::cout<<"DONE"<<std::endl;
			std::cout<<"_____________________________________________________"<<std::endl;
			}


			/*********************************************************************************************
			 * OUTPUT
			 *********************************************************************************************/
			/* At each time step: **Statistics**
			 * t_new,
			 * dt,
			 * fixed point iterations,
			 * newton iterations per fixed point iteration,
			 * total newton terations
			 */
			if(helper.rank()==0){
			std::string statistics_file = pathName;
			statistics_file +=fileName;
			statistics_file +="_statistics";
			statistics_file += ".txt";
			property.ReportStatistics( 	statistics_file,
										time*CharacteristicValues::t_c,
										dt*CharacteristicValues::t_c,
										newton_iterations,
										clock_time_elapsed );
			}

			/* At t_OP
			 *
			 */
			if( ( time+dt > t_OP*opcount - dt_min ) and ( time+dt <= t_OP*opcount )  )
			{
				// POST PROCESS FOR NEW OUTPUT
				postprocess.evaluate();
#ifdef PLOT_VELOCITIES
				slpproj.apply();
				gasvelocity.evaluate();
				watervelocity.evaluate();
#endif

				// WRITE OUTPUT
				vtkSequenceWriter.write(time,Dune::VTK::appendedraw);

				if(helper.rank()==0){
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< " OUTPUT WRITTEN " << opcount << std::endl;
				std::cout<< " ******************************************************************* " << std::endl;
				std::cout<< std::flush;
				}

				timecount = time;
				opcount = opcount+1;
			}

			//	PREPARE FOR NEXT TIME INTEGRATION
			//	1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
			u_old = u_new;
			//	2. ADVANCE TIME:
			time += dt;

			if(helper.rank()==0){
			std::cout<<" "<< std::endl;
			std::cout<< " time = " << time*property.characteristicValue.t_c ;
			std::cout<< std::flush;
			}

			if( adaptive_time_control ){
				if(newton_iterations>maxAllowableIterations){
					dt=std::max(dt*0.9,dt_min);
				}
				else if(newton_iterations<=minAllowableIterations){
					dt=std::min(dt*1.1,dt_max);
				}
			}
			else{
				dt = dtstart;
			}

			if(helper.rank()==0){
			std::cout << " , time+dt = " << (time + dt)*property.characteristicValue.t_c
					  << " , opTime = "  << t_OP * opcount * property.characteristicValue.t_c ;
			std::cout<< std::flush;
			}

			if( time + dt  > t_OP * opcount){
				dtLast = dt;
				dt 	 = t_OP * opcount - time ;

				if(helper.rank()==0){
				std::cout<< " , because timeNext > opNext , dt set to : " << dt*property.characteristicValue.t_c << std::endl;
				std::cout<< std::flush;
				}

				dtFlag = 0;
			}
			dtFlag += 1;

			if(helper.rank()==0){
			std::cout<< " , dt  : " << dt*property.characteristicValue.t_c << std::endl;
			std::cout<<" "<< std::endl;
			std::cout << " READY FOR NEXT ITERATION. " << std::endl;
			std::cout<< std::flush;
			}
		}

}

