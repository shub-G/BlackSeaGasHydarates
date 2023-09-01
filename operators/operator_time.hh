template < class GV, class Params >
class TimeOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume< TimeOperator<GV,Params> >,
    public Dune::PDELab::NumericalJacobianVolume	 < TimeOperator<GV,Params> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV& gv;
	const Params&	  param;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double pi 		= atan(1.)*4 ;
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

public:
	  // pattern assembly flags
	  enum { doPatternVolume = true };

	  // residual assembly flags
	  enum { doAlphaVolume = true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor remembers parameters
	  TimeOperator( const GV& gv_, const Params& param_ )
	  :  gv(gv_), param(param_)
	  {
		  Xc_K 		= param.characteristicValue.permeability_c;
		  Xc_mu 	= param.characteristicValue.viscosity_c;
		  Xc_rho 	= param.characteristicValue.density_c;
		  Xc_kth 	= param.characteristicValue.thermalconductivity_c;
		  Xc_C 		= param.characteristicValue.specificheat_c;
		  Xc_D 		= param.characteristicValue.dispersivity_c;
		  Xc_P 		= param.characteristicValue.P_c;
		  Xc_T		= param.characteristicValue.T_c;
		  Xc_t 		= param.characteristicValue.t_c;
		  Xc_x 		= param.characteristicValue.x_c;
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)const{

	        typedef typename LFSU::template Child< Indices::PVId_Pw	  >::Type LFS_Pw ;
	        const LFS_Pw&	lfs_Pw	= lfsu.template child< Indices::PVId_Pw  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sg   >::Type LFS_Sg;
	        const LFS_Sg&	lfs_Sg	= lfsu.template child< Indices::PVId_Sg  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_XCH4 >::Type LFS_XCH4;
	        const LFS_XCH4& lfs_XCH4= lfsu.template child< Indices::PVId_XCH4>() ;
	        typedef typename LFSU::template Child< Indices::PVId_YH2O >::Type LFS_YH2O;
	        const LFS_YH2O& lfs_YH2O= lfsu.template child< Indices::PVId_YH2O>() ;
	        typedef typename LFSU::template Child< Indices::PVId_Xc   >::Type LFS_Xc ;
	        const LFS_Xc& 	lfs_Xc	= lfsu.template child< Indices::PVId_Xc  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sh   >::Type LFS_Sh ;
	        const LFS_Sh&	lfs_Sh	= lfsu.template child< Indices::PVId_Sh  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_T    >::Type LFS_T  ;
	        const LFS_T& 	lfs_T 	= lfsu.template child< Indices::PVId_T   >() ;

			// Reference to cell
	        const auto& cell = eg.entity();
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);

	        // get geometry
	        auto geo = eg.geometry();

			// dimension
			const auto dim = geo.mydimension;

	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_center_global = geo.global(cell_center_local);
	        auto cell_volume = geo.volume();
	        cell_volume *= param.mesh.volumeFactor(cell_center_global[0]);

			// compute PVs at local center
			double Pw   = x(lfs_Pw  ,0);
			double Sg   = x(lfs_Sg  ,0);
			double XCH4 = x(lfs_XCH4,0);
			double YH2O = x(lfs_YH2O,0);
			double Xc   = x(lfs_Xc  ,0);
			double Sh   = x(lfs_Sh  ,0);
			double T    = x(lfs_T   ,0);

			auto porosity = param.soil.SedimentPorosity( cell,cell_center_local );
#ifdef COMPACTION
			porosity *= param.parameter.CompactionFunction( cell_center_global );
#endif
			auto permeability = param.soil.SedimentPermeability( cell,cell_center_local )
							  * param.hydraulicProperty.PermeabilityScalingFactor( cell,cell_center_local, Sh, porosity );
			auto S = Xc * (param.salt.MolarMass()/param.water.MolarMass());
			auto Sw = 1.-Sh-Sg;
			auto Pc = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local, Sw, Sh, porosity );
			auto Pg = Pw + Pc;
			auto Peff = ( Pg*Sg + Pw*Sw )/(1.-Sh);
			auto zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );

			// compute properties at local center
			// density
			auto rho_mass_g = param.gas.Density( T*Xc_T, Pg*Xc_P, zCH4 );
			auto rho_mass_w = param.water.Density( T*Xc_T, Pw*Xc_P, S );
			auto rho_mass_h  = param.hydrate.Density();
			auto rho_mass_s  = param.soil.Density();
			// specific heat capacities
			auto Cv_g = param.gas.Cv( T*Xc_T, Pg*Xc_P, zCH4 );
			auto Cv_w = param.water.Cv( T*Xc_T, Pw*Xc_P, S );
			auto rhoCv_g = rho_mass_g * Cv_g;
			auto rhoCv_w = rho_mass_w * Cv_w;
			auto rhoCv_h = rho_mass_h * param.hydrate.Cv( T*Xc_T, Peff*Xc_P );
			auto rhoCv_s = rho_mass_s * param.soil.Cv();
			auto Cv_eff = (1.-porosity) * rhoCv_s
						+ porosity * Sg * rhoCv_g
						+ porosity * Sw * rhoCv_w
						+ porosity * Sh * rhoCv_h;


			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_g * (1.-YH2O) * Sg ;
			tmp += porosity * rho_mass_w * XCH4 * Sw ;
//			std::cout<< "alpha-time 0: " << tmp << std::endl;
			r.accumulate(lfs_Pw, 0, +tmp*cell_volume);

			// H2O-component-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_g * YH2O * Sg ;
			tmp += porosity * rho_mass_w * (1.-Xc-XCH4) * Sw ;
			r.accumulate(lfs_Sg, 0, +tmp*cell_volume);

			// SALT-component-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_w * Xc * Sw ;
			r.accumulate(lfs_Xc, 0, +tmp*cell_volume);

			// HYD-phase-wise mass-balance
			tmp = 0.;
			tmp += porosity * rho_mass_h * Sh;
			r.accumulate(lfs_Sh, 0, +tmp*cell_volume);

			// ENERGY balance
			auto T_ref = param.parameter.ReferenceTemperature();
			tmp = 0.;
			tmp += Cv_eff * (T-T_ref);
			r.accumulate(lfs_T , 0, +tmp*cell_volume);

	  }

};
