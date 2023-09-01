/*********************************************************
 * EVALUATE OUTPUT VARIABLES
 *********************************************************/

template< class GV,
		  class Params,
		  class Evaluation_Pw,
		  class Evaluation_Sg,
		  class Evaluation_XCH4,
		  class Evaluation_YH2O,
		  class Evaluation_Xc,
		  class Evaluation_Sh,
		  class Evaluation_T,
		  class GFS_PP, typename U_pp >
class PostProcess{
private:
	const GV&  gv;
	const Params& param;
	Evaluation_Pw    *evaluation_Pw;
	Evaluation_Sg    *evaluation_Sg;
	Evaluation_XCH4  *evaluation_XCH4;
	Evaluation_YH2O  *evaluation_YH2O;
	Evaluation_Xc    *evaluation_Xc;
	Evaluation_Sh    *evaluation_Sh;
	Evaluation_T     *evaluation_T ;
	GFS_PP gfs_pp;
	U_pp *u_pp;

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

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PostProcess(	const 	GV& gv_,
					const Params& 	param_,
					Evaluation_Pw 	*evaluation_Pw_,
					Evaluation_Sg 	*evaluation_Sg_,
					Evaluation_XCH4 *evaluation_XCH4_,
					Evaluation_YH2O *evaluation_YH2O_,
					Evaluation_Xc 	*evaluation_Xc_,
					Evaluation_Sh 	*evaluation_Sh_,
					Evaluation_T  	*evaluation_T_ ,
					GFS_PP	gfs_pp_,
					U_pp	*u_pp_)
	: gv(gv_),
	  param(param_),
	  evaluation_Pw(evaluation_Pw_),
	  evaluation_Sg(evaluation_Sg_),
	  evaluation_XCH4(evaluation_XCH4_),
	  evaluation_YH2O(evaluation_YH2O_),
	  evaluation_Xc(evaluation_Xc_),
	  evaluation_Sh(evaluation_Sh_),
	  evaluation_T (evaluation_T_ ),
	  gfs_pp(gfs_pp_),
	  u_pp(u_pp_)
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
	  }

	virtual ~PostProcess()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS_PP > LFS_PP;
		LFS_PP lfs_pp(gfs_pp);
		typedef typename LFS_PP::template Child<Indices::SVId_Pg>::Type LFS_PP_Pg;
		const LFS_PP_Pg& lfs_pp_Pg = lfs_pp.template child<Indices::SVId_Pg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pw>::Type LFS_PP_Pw;
		const LFS_PP_Pw& lfs_pp_Pw = lfs_pp.template child<Indices::SVId_Pw>();
		typedef typename LFS_PP::template Child<Indices::SVId_Pc>::Type LFS_PP_Pc;
		const LFS_PP_Pc& lfs_pp_Pc = lfs_pp.template child<Indices::SVId_Pc>();

		typedef typename LFS_PP::template Child<Indices::SVId_Sg>::Type LFS_PP_Sg;
		const LFS_PP_Sg& lfs_pp_Sg = lfs_pp.template child<Indices::SVId_Sg>();
		typedef typename LFS_PP::template Child<Indices::SVId_Sw>::Type LFS_PP_Sw;
		const LFS_PP_Sw& lfs_pp_Sw = lfs_pp.template child<Indices::SVId_Sw>();
		typedef typename LFS_PP::template Child<Indices::SVId_Sh>::Type LFS_PP_Sh;
		const LFS_PP_Sh& lfs_pp_Sh = lfs_pp.template child<Indices::SVId_Sh>();

		typedef typename LFS_PP::template Child<Indices::SVId_XCH4>::Type LFS_PP_XCH4;
		const LFS_PP_XCH4& lfs_pp_XCH4 = lfs_pp.template child<Indices::SVId_XCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_XH2O>::Type LFS_PP_XH2O;
		const LFS_PP_XH2O& lfs_pp_XH2O = lfs_pp.template child<Indices::SVId_XH2O>();
		typedef typename LFS_PP::template Child<Indices::SVId_YCH4>::Type LFS_PP_YCH4;
		const LFS_PP_YCH4& lfs_pp_YCH4 = lfs_pp.template child<Indices::SVId_YCH4>();
		typedef typename LFS_PP::template Child<Indices::SVId_YH2O>::Type LFS_PP_YH2O;
		const LFS_PP_YH2O& lfs_pp_YH2O = lfs_pp.template child<Indices::SVId_YH2O>();

		typedef typename LFS_PP::template Child<Indices::SVId_Xc>::Type LFS_PP_Xc;
		const LFS_PP_Xc& lfs_pp_Xc = lfs_pp.template child<Indices::SVId_Xc>();

		typedef typename LFS_PP::template Child<Indices::SVId_T>::Type LFS_PP_T;
		const LFS_PP_T& lfs_pp_T = lfs_pp.template child<Indices::SVId_T>();

		typedef typename LFS_PP::template Child<Indices::SVId_K>::Type LFS_PP_K;
		const LFS_PP_K& lfs_pp_K = lfs_pp.template child<Indices::SVId_K>();

		typedef typename LFS_PP::template Child<Indices::SVId_z>::Type LFS_PP_zCH4;
		const LFS_PP_zCH4& lfs_pp_zCH4 = lfs_pp.template child<Indices::SVId_z>();

		typedef typename LFS_PP::template Child<Indices::SVId_Peq>::Type LFS_PP_Peq;
		const LFS_PP_Peq& lfs_pp_Peq = lfs_pp.template child<Indices::SVId_Peq>();

		typedef typename LFS_PP::template Child<Indices::SVId_por>::Type LFS_PP_por;
		const LFS_PP_por& lfs_pp_por = lfs_pp.template child<Indices::SVId_por>();

		typedef Dune::PDELab::LFSIndexCache<LFS_PP> LFSCache_PP;
		LFSCache_PP lfs_cache_pp(lfs_pp);
		typedef typename U_pp::template LocalView<LFSCache_PP> VectorView_PP;
		VectorView_PP u_pp_view( (*u_pp) );

		typedef typename LFS_PP::template Child<0>::Type::Traits::FiniteElementType::
		  Traits::LocalBasisType::Traits::RangeFieldType RF;

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

	        RF Pw=0.;
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
	        RF T=0.;
	        evaluation_T->evalFunction(cell,cell_center_local,&T);

	        lfs_pp.bind(*self);
			lfs_cache_pp.update();
			u_pp_view.bind(lfs_cache_pp);
	        std::vector<double> ul_pp(lfs_pp.size());
	        for(int i = 0. ; i < lfs_pp.size() ; i++){
	        	ul_pp[i] = 0.;
	        }

			auto porosity = param.soil.SedimentPorosity( cell,cell_center_local );
#ifdef COMPACTION
			porosity *= param.parameter.CompactionFunction( cell_center_global );
#endif
			auto K = param.soil.SedimentPermeability( cell,cell_center_local )
				   * param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, Sh, porosity );

			RF S = Xc*(param.salt.MolarMass()/param.water.MolarMass());
			RF Sw = 1.-Sh-Sg;
			RF Pc = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, Sh, porosity );
			RF Pg = Pw + Pc;
			RF zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );
			RF XH2O = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,Xc);
			RF YCH4 = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,Xc,zCH4);
			RF Peq = param.kinetics.EquilibriumPressure(T*Xc_T,S);


			ul_pp[lfs_pp_Pg.localIndex(0)] = Pg*Xc_P ;
	        ul_pp[lfs_pp_Pw.localIndex(0)] = Pw*Xc_P ;
	        ul_pp[lfs_pp_Pc.localIndex(0)] = Pc*Xc_P ;
	        ul_pp[lfs_pp_Sg.localIndex(0)] = Sg ;
	        ul_pp[lfs_pp_Sw.localIndex(0)] = Sw ;
	        ul_pp[lfs_pp_Sh.localIndex(0)] = Sh ;
	        ul_pp[lfs_pp_XCH4.localIndex(0)] = XCH4 ;
	        ul_pp[lfs_pp_XH2O.localIndex(0)] = XH2O ;
	        ul_pp[lfs_pp_YCH4.localIndex(0)] = YCH4 ;
	        ul_pp[lfs_pp_YH2O.localIndex(0)] = YH2O ;
	        ul_pp[lfs_pp_Xc.localIndex(0)] = Xc ;
	        ul_pp[lfs_pp_T.localIndex(0)] = T*Xc_T ;
	        ul_pp[lfs_pp_K.localIndex(0)] = K*Xc_K ;
	        ul_pp[lfs_pp_zCH4.localIndex(0)] = zCH4 ;
	        ul_pp[lfs_pp_Peq.localIndex(0)] = Peq*Xc_P ;
	        ul_pp[lfs_pp_por.localIndex(0)] = porosity ;

			u_pp_view.write( ul_pp );
			u_pp_view.commit();
			u_pp_view.unbind();

		}//END:iterate over each volume

	}
};
