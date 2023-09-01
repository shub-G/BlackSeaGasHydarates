template< class GV,
		  class Params,
		  class Evaluation_PROJP,
		  class Evaluation_Pw,
		  class Evaluation_Sg,
		  class Evaluation_XCH4,
		  class Evaluation_YH2O,
		  class Evaluation_Xc,
		  class Evaluation_Sh,
		  class Evaluation_T,
		  class GFS_V, typename U_V >
class PhaseVelocity{
private:
	const GV&  gv;
	const Params& param;
	double *time;
	double *dt;
	const int phase_tag;
	Evaluation_PROJP *evaluation_P;
	Evaluation_Pw    *evaluation_Pw;
	Evaluation_Sg    *evaluation_Sg;
	Evaluation_XCH4  *evaluation_XCH4;
	Evaluation_YH2O  *evaluation_YH2O;
	Evaluation_Xc    *evaluation_Xc;
	Evaluation_Sh    *evaluation_Sh;
	Evaluation_T     *evaluation_T ;
	GFS_V gfs_v;
	U_V *u_v;

	double Xc_K ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_kth;
	double Xc_C;
	double Xc_D;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;
	double Xc_grav;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PhaseVelocity(	const 	GV& gv_,
					const Params& param_,
					Evaluation_PROJP *evaluation_P_,
					Evaluation_Pw 	 *evaluation_Pw_,
					Evaluation_Sg 	 *evaluation_Sg_,
					Evaluation_XCH4  *evaluation_XCH4_,
					Evaluation_YH2O  *evaluation_YH2O_,
					Evaluation_Xc 	 *evaluation_Xc_,
					Evaluation_Sh 	 *evaluation_Sh_,
					Evaluation_T  	 *evaluation_T_ ,
					GFS_V	gfs_v_,
					const int phase_tag_,
					U_V	*u_v_,
					double *time_,
					double *dt_ )
	: gv(gv_),
	  param(param_),
	  evaluation_P(evaluation_P_),
	  evaluation_Pw(evaluation_Pw_),
	  evaluation_Sg(evaluation_Sg_),
	  evaluation_XCH4(evaluation_XCH4_),
	  evaluation_YH2O(evaluation_YH2O_),
	  evaluation_Xc(evaluation_Xc_),
	  evaluation_Sh(evaluation_Sh_),
	  evaluation_T (evaluation_T_ ),
	  gfs_v(gfs_v_),
	  phase_tag(phase_tag_),
	  u_v(u_v_),
	  time(time_),
	  dt(dt_)
	{
		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_kth 	= param.characteristicValue.thermalconductivity_c;
		  Xc_C 		= param.characteristicValue.specificheat_c;
		  Xc_D		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T 		= param.characteristicValue.T_c;
		  Xc_t 		= param.characteristicValue.t_c;
		  Xc_x 		= param.characteristicValue.x_c;
		  Xc_grav 	= param.characteristicValue.X_gravity;
	  }

	virtual ~PhaseVelocity()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS_V > LFS_V;
		LFS_V lfs_v(gfs_v);
		typedef Dune::PDELab::LFSIndexCache<LFS_V> LFSCache_V;
		LFSCache_V lfs_cache_v(lfs_v);
		typedef typename U_V::template LocalView<LFSCache_V> VectorView_V;
		VectorView_V u_v_view( (*u_v) );

		using RF = typename LFS_V::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
		using RangeType = typename LFS_V::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
		using JacobianType = typename LFS_V::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
		typedef typename LeafIterator::Entity::Geometry EG;

		// Loop over each volume
		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_center_global = geo.global(cell_center_local);
	        auto cell_volume = geo.volume();

	        Dune::FieldVector<RF,dim> gradP(0.);
	        evaluation_P->evalGradient(cell,cell_center_local,&gradP);
	        RF P_proj/*ndim*/=0.;
	        evaluation_P->evalFunction(cell,cell_center_local,&P_proj);
	        RF Pw/*ndim*/=0.;
	        evaluation_Pw->evalFunction(cell,cell_center_local,&Pw);
	        RF Sg=0.;
	        evaluation_Sg->evalFunction(cell,cell_center_local,&Sg);
	        RF XCH4=0.;
	        evaluation_XCH4->evalFunction(cell,cell_center_local,&XCH4);
	        RF YH2O=0.;
	        evaluation_YH2O->evalFunction(cell,cell_center_local,&YH2O);
	        RF Xc=0.;
	        evaluation_Xc->evalFunction(cell,cell_center_local,&Xc);
	        RF Sh=0.;
	        evaluation_Sh->evalFunction(cell,cell_center_local,&Sh);
	        RF T/*ndim*/=0.;
	        evaluation_T->evalFunction(cell,cell_center_local,&T);


	        lfs_v.bind(*self);
			lfs_cache_v.update();
			u_v_view.bind(lfs_cache_v);
	        std::vector<double> ul_v(lfs_v.size());
	        for(int i = 0. ; i < lfs_v.size() ; i++){
	        	ul_v[i] = 0.;
	        }


			auto porosity = param.soil.SedimentPorosity( cell,cell_center_local );
#ifdef COMPACTION
			porosity *= param.parameter.CompactionFunction( cell_center_global );
#endif
			auto K/*ndim*/ = param.soil.SedimentPermeabilityVector( cell,cell_center_local );
			K *= param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, Sh, porosity );

			auto S = Xc * (param.salt.MolarMass()/param.water.MolarMass());
			auto Sw = 1.-Sh-Sg;
			auto Pc/*ndim*/ = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, Sh, porosity );
			auto Pg/*ndim*/ = Pw + Pc;
			auto zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );

			RF phase_flag = 0.;
			RF lambda = 1.;
			RF rho/*ndim*/ = 1.;
			if(phase_tag==Indices::phaseId_g){
				auto mu = param.gas.DynamicViscosity( T*Xc_T, Pg*Xc_P );
				auto kr = param.hydraulicProperty.krg( cell, cell_center_local, Sw, Sh );
				lambda = kr/mu;
				rho    = param.gas.Density( T*Xc_T, Pg*Xc_P, zCH4 );
				if(Sg>0.){
					phase_flag = 1.;
				}
			}else if(phase_tag==Indices::phaseId_w){
				auto mu = param.water.DynamicViscosity( T*Xc_T, Pw*Xc_P, S );
				auto kr = param.hydraulicProperty.krw( cell, cell_center_local, Sw, Sh );
				lambda = kr/mu;
				rho    = param.water.Density( T*Xc_T, Pw*Xc_P, S );
				if(Sw>0.){
					phase_flag = 1.;
				}
			}else{
				std::cout<< "WARNING: Invalid phase ID in: " << __FILE__ << std::endl;
				exit(0);
			}

			auto gravity/*m/s^2*/ = param.parameter.g() ;

			for(int d=0; d<dim; d++){
				ul_v[lfs_v.child(d).localIndex(0)] = - phase_flag * K[d]*Xc_K * (lambda/Xc_mu) * ( gradP[d]*Xc_P/Xc_x + (rho*Xc_rho)*gravity[d] );
//				ul_v[lfs_v.child(d).localIndex(0)] = - K[d]*Xc_K * (lambda/Xc_mu) * ( gradP[d]*Xc_P/Xc_x + (rho*Xc_rho)*gravity[d] );
			}

			u_v_view.write( ul_v );
			u_v_view.commit();
			u_v_view.unbind();

		}//END:iterate over each volume

	}
};

