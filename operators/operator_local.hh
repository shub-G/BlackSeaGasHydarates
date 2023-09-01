template <	class GV, class Params, class BC >
class LocalOperator :
  public Dune::PDELab::NumericalJacobianApplyVolume		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianVolume			< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianApplySkeleton	< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianSkeleton		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianApplyBoundary	< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::NumericalJacobianBoundary		< LocalOperator<GV,Params,BC> >,
  public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags,
  public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
    const GV& gv;
	const Params&	  param;
	const BC&	 	  bc;
	constexpr static double eps 	= 1.0e-6;
	constexpr static double pi 		= atan(1.)*4 ;
	double Xc_conv_m ;
	double Xc_conv_h ;
	double Xc_source_m ;
	double Xc_source_h ;
	double Xc_diff_m ;
	double Xc_diff_h ;
	double Xc_grav ;
	double Xc_vs ;
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
	  enum { doPatternVolume	= true };
	  enum { doPatternSkeleton	= true };

	  // residual assembly flags
	  enum { doAlphaVolume  	= true };
	  enum { doAlphaSkeleton	= true };
	  enum { doAlphaBoundary	= true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor stores parameters
	  LocalOperator( const GV& 		 gv_,
				  	 const Params&	 param_,
					 const BC& 	 	 bc_	)
	  :  gv(gv_),
		 param( param_ ),
		 bc( bc_ )
	  {
		  Xc_conv_m 	= param.characteristicValue.X_convective_mass;
		  Xc_conv_h 	= param.characteristicValue.X_convective_heat;
		  Xc_source_m 	= param.characteristicValue.X_source_mass;
		  Xc_source_h 	= param.characteristicValue.X_source_heat;
		  Xc_diff_m 	= param.characteristicValue.X_diffusive_mass;
		  Xc_diff_h 	= param.characteristicValue.X_diffusive_heat;
		  Xc_grav 		= param.characteristicValue.X_gravity;
		  Xc_vs			= param.characteristicValue.X_solidvelocity;

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

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
	  {
	        typedef typename LFSU::template Child< Indices::PVId_Pw	 >::Type LFS_Pw ;
	        const LFS_Pw&	lfs_Pw	= lfsu.template child< Indices::PVId_Pw  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sg	 >::Type LFS_Sg;
	        const LFS_Sg& 	lfs_Sg	= lfsu.template child< Indices::PVId_Sg  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_XCH4>::Type LFS_XCH4;
	        const LFS_XCH4&	lfs_XCH4= lfsu.template child< Indices::PVId_XCH4>() ;
	        typedef typename LFSU::template Child< Indices::PVId_YH2O>::Type LFS_YH2O;
	        const LFS_YH2O&	lfs_YH2O= lfsu.template child< Indices::PVId_YH2O>() ;
	        typedef typename LFSU::template Child< Indices::PVId_Xc  >::Type LFS_Xc ;
	        const LFS_Xc&	lfs_Xc	= lfsu.template child< Indices::PVId_Xc  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_Sh  >::Type LFS_Sh ;
	        const LFS_Sh&	lfs_Sh  = lfsu.template child< Indices::PVId_Sh  >() ;
	        typedef typename LFSU::template Child< Indices::PVId_T   >::Type LFS_T  ;
	        const LFS_T&	lfs_T   = lfsu.template child< Indices::PVId_T   >() ;

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
	        auto cell_center_global = cell.geometry().global(cell_center_local);
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
			auto Peff = (Pg*Sg+Pw*Sw)/(1.-Sh);
			auto zCH4 = param.eos.EvaluateCompressibilityFactor( T*Xc_T,Pg*Xc_P );

			// compute source terms
			auto g_CH4  = param.kinetics.GasGenerationRate( T*Xc_T,
														    Pg*Xc_P,
														    Sh,
														    Sw,
														    XCH4,
														    zCH4,
														    S,
														    porosity,
														    permeability*Xc_K);
			auto g_H2O  = param.kinetics.WaterGenerationRate( g_CH4 );
			auto g_HYD  = param.kinetics.HydrateDissociationRate( g_CH4 );
			auto g_salt = param.salt.Source();
			auto Q_diss = param.kinetics.HeatOfDissociation( g_CH4, T*Xc_T );



			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp -= Xc_source_m * g_CH4 ;
//			std::cout<< "alpha_vol 0: " << tmp << std::endl;
			r.accumulate(lfs_Pw, 0, +tmp*cell_volume);

			// H2O-component-wise mass-balance
			tmp = 0.;
			tmp -= Xc_source_m * g_H2O ;
//			std::cout<< "alpha_vol 1: " << tmp << std::endl;
			r.accumulate(lfs_Sg, 0, +tmp*cell_volume);

			// SALT-component-wise mass-balance
			tmp = 0.;
			tmp -= Xc_source_m * g_salt ;
//			std::cout<< "alpha_vol 2: " << tmp << std::endl;
			r.accumulate(lfs_Xc, 0, +tmp*cell_volume);

			// HYD-phase-wise mass-balance
			tmp = 0.;
			tmp -= Xc_source_m * g_HYD ;
//			std::cout<< "alpha_vol 3: " << tmp << std::endl;
			r.accumulate(lfs_Sh, 0., +tmp*cell_volume);

			// ENERGY balance
			tmp = 0.;
			tmp -= Xc_source_h * Q_diss;
//			std::cout<< "alpha_vol 4: " << tmp << std::endl;
			r.accumulate(lfs_T , 0., +tmp*cell_volume);


//			double porosity0 = param.soil.SedimentPorosity(cell,cell_center_local);
//			double Pw_PSF0 = param.parameter.PSF_Pressure()/param.characteristicValue.P_c; /*ndim*/
//			double Pw0 = Pw_PSF0 + 1000.*10.*(param.mesh.z_PSFC(cell_center_global[0])-cell_center_global[1])
//					   * param.characteristicValue.x_c/param.characteristicValue.P_c; /*ndim*/
//			double Pg0 = Pw0 + Pc; /*ndim*/
//
//			double T_PSF0 	= param.parameter.PSF_Temperature()/param.characteristicValue.T_c; /*ndim*/
//			double T_grad0	= param.parameter.RegionalThermalGradient()
//							* param.characteristicValue.x_c/param.characteristicValue.T_c; /*ndim*/
//			double T0 = T_PSF0;
//			if(cell_center_global[1]<param.mesh.z_PSFC(cell_center_global[0])+eps){
//				T0 +=  T_grad0*(param.mesh.z_PSFC(cell_center_global[0])-cell_center_global[1]);
//			}
			// NCP -> water phase
			tmp = 0.;
//			auto XH2O_alg = param.mixture.XH2O(YH2O,T*Xc_T,Pg*Xc_P,Xc);
//			if( ( Sw - ( 1. - XCH4 - XH2O_alg - Xc ) ) > 0. ){//active set.
//				tmp += 1. - XCH4 - XH2O_alg - Xc;//Active => phase is present => summation condition holds
////				std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
//			}
			auto YH2O_eqb = param.water.EquilibriumYH2O(T*Xc_T,Pg*Xc_P,S);
			if( ( Sw - ( YH2O_eqb - YH2O ) ) > 0. ){//active set.
				tmp +=  YH2O_eqb - YH2O ;//Active => phase is present => summation condition holds
//				std::cout<< "alpha_vol XCH4: " << tmp << std::endl;
			}
			else{
				tmp += Sw; // inactive set. Inactive => phase is absent => Sw=0
			}
			r.accumulate(lfs_XCH4 , 0., +tmp*cell_volume);

			// NCP -> gas phase
			tmp = 0.;
//			auto YCH4_alg = param.mixture.YCH4(XCH4,T*Xc_T,Pg*Xc_P,Xc,zCH4);
//			if( ( Sg - ( 1. - YCH4_alg - YH2O ) ) > 0. ){//active set.
			auto XCH4_eqb = param.gas.EquilibriumXCH4(T*Xc_T,Pg*Xc_P,zCH4,S);
			if( ( Sg - ( XCH4_eqb - XCH4 ) ) > 0. ){//active set.
				tmp += XCH4_eqb - XCH4;//1. - YCH4_alg - YH2O ;//Active => phase is present => summation condition holds
//				std::cout<< "alpha_vol YH2O: " << tmp << std::endl;
			}else{
				tmp += Sg;// inactive set. Inactive => phase is absent => Sg=0
			}
			r.accumulate(lfs_YH2O , 0., +tmp*cell_volume);

	  }

	  // skeleton integral depending on test and ansatz functions
	  // each face is only visited ONCE!
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_skeleton (const IG& ig,
			  	  	  	   const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
						   const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
						   R& r_s, R& r_n) const
	  {
			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_s;
			const LFS_Pw_s& lfs_Pw_s = lfsu_s.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_s;
			const LFS_Sg_s& lfs_Sg_s = lfsu_s.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_s;
			const LFS_XCH4_s& lfs_XCH4_s = lfsu_s.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_s;
			const LFS_YH2O_s& lfs_YH2O_s = lfsu_s.template child<Indices::PVId_YH2O>();
			typedef typename LFSU::template Child<Indices::PVId_Xc>::Type LFS_Xc_s;
			const LFS_Xc_s& lfs_Xc_s = lfsu_s.template child<Indices::PVId_Xc>();
			typedef typename LFSU::template Child<Indices::PVId_Sh>::Type LFS_Sh_s;
			const LFS_Sh_s& lfs_Sh_s = lfsu_s.template child<Indices::PVId_Sh>();
			typedef typename LFSU::template Child<Indices::PVId_T >::Type LFS_T_s ;
			const LFS_T_s&  lfs_T_s  = lfsu_s.template child<Indices::PVId_T >();

			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_n;
			const LFS_Pw_n& lfs_Pw_n = lfsu_n.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_n;
			const LFS_Sg_n& lfs_Sg_n = lfsu_n.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_n;
			const LFS_XCH4_n& lfs_XCH4_n = lfsu_n.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_n;
			const LFS_YH2O_n& lfs_YH2O_n = lfsu_n.template child<Indices::PVId_YH2O>();
			typedef typename LFSU::template Child<Indices::PVId_Xc>::Type LFS_Xc_n;
			const LFS_Xc_n& lfs_Xc_n = lfsu_n.template child<Indices::PVId_Xc>();
			typedef typename LFSU::template Child<Indices::PVId_Sh>::Type LFS_Sh_n;
			const LFS_Sh_n& lfs_Sh_n = lfsu_n.template child<Indices::PVId_Sh>();
			typedef typename LFSU::template Child<Indices::PVId_T >::Type LFS_T_n ;
			const LFS_T_n&  lfs_T_n  = lfsu_n.template child<Indices::PVId_T >();

	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
	        const auto& cell_outside = ig.outside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);
			int outside_cell_number = indexSet.index(cell_outside);

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = geo.mydimension;
	        auto geo_inside  = cell_inside.geometry();
	        auto geo_outside = cell_outside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto ref_el_outside = referenceElement(geo_outside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto outside_cell_center_local 	= ref_el_outside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();
	        auto outside_cell_center_global = geo_outside.center();

	        // distance of cell centers
	        auto d = outside_cell_center_global;
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);
	        face_volume *= param.mesh.areaFactor( inside_cell_center_global[0],
	        									  face_center_global[0],
												  normal);

			// compute primary vars at local self and neighbour centers
	        double Pw_s   = x_s(lfs_Pw_s,  0);
			double Sg_s   = x_s(lfs_Sg_s,  0);
			double XCH4_s = x_s(lfs_XCH4_s,0);
			double YH2O_s = x_s(lfs_YH2O_s,0);
			double Xc_s   = x_s(lfs_Xc_s,  0);
			double Sh_s   = x_s(lfs_Sh_s,  0);
			double T_s    = x_s(lfs_T_s ,  0);

	        double Pw_n   = x_n(lfs_Pw_n,  0);
			double Sg_n   = x_n(lfs_Sg_n,  0);
			double XCH4_n = x_n(lfs_XCH4_n,0);
			double YH2O_n = x_s(lfs_YH2O_n,0);
			double Xc_n   = x_n(lfs_Xc_n,  0);
			double Sh_n   = x_n(lfs_Sh_n,  0);
			double T_n    = x_n(lfs_T_n ,  0);

			// compute properties at local cell center and neighbour cell center

			// Hydraulic properties for self and neighbour
			// 1. porosity
			auto porosity_s = param.soil.SedimentPorosity( cell_inside,inside_cell_center_local );
			auto porosity_n = param.soil.SedimentPorosity( cell_outside,outside_cell_center_local );
#ifdef COMPACTION
			porosity_s *= param.parameter.CompactionFunction( inside_cell_center_global );
			porosity_n *= param.parameter.CompactionFunction( outside_cell_center_global );
#endif
			// 2. abs. permeability
			auto permeability_s =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local) ;
			permeability_s *= param.hydraulicProperty.PermeabilityScalingFactor( cell_inside,inside_cell_center_local, Sh_s, porosity_s );
			permeability_s = std::abs(permeability_s);
			auto permeability_n =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local);
			permeability_n *= param.hydraulicProperty.PermeabilityScalingFactor( cell_outside,outside_cell_center_local, Sh_s, porosity_s );
			permeability_n = std::abs(permeability_n);
			// 3. tortuosity
			auto tau_s = param.soil.Tortuosity( porosity_s*(1.-Sh_s) );
			auto tau_n = param.soil.Tortuosity( porosity_n*(1.-Sh_n) );

			// 4. secondary variables
			auto S_s = Xc_s * (param.salt.MolarMass()/param.water.MolarMass());
			auto S_n = Xc_n * (param.salt.MolarMass()/param.water.MolarMass());
			auto Sw_s = 1.-Sh_s-Sg_s;
			auto Sw_n = 1.-Sh_n-Sg_n;
			auto Pc_s = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_s, Sh_s, porosity_s );
			auto Pc_n = param.hydraulicProperty.CapillaryPressure( cell_outside,outside_cell_center_local, Sw_n, Sh_n, porosity_n );
			auto Pg_s = Pw_s + Pc_s;
			auto Pg_n = Pw_n + Pc_n;
			auto Peff_s = ( Pg_s*Sg_s + Pw_s*Sw_s )/( 1.-Sh_s );
			auto Peff_n = ( Pg_n*Sg_n + Pw_n*Sw_n )/( 1.-Sh_n );
			auto zCH4_s = param.eos.EvaluateCompressibilityFactor( T_s*Xc_T,Pg_s*Xc_P );
			auto zCH4_n = param.eos.EvaluateCompressibilityFactor( T_n*Xc_T,Pg_n*Xc_P );

			// 5. relative permeability
			// self
			auto kr_g_s = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_s, Sh_s );
			auto kr_w_s = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_s, Sh_s );
			// neighbour
			auto kr_g_n = param.hydraulicProperty.krg( cell_outside, outside_cell_center_local, Sw_n, Sh_n );
			auto kr_w_n = param.hydraulicProperty.krw( cell_outside, outside_cell_center_local, Sw_n, Sh_n );

			// 6. dynamic viscosity
			// self
			auto mu_g_s = param.gas.DynamicViscosity( T_s*Xc_T, Pg_s*Xc_P );
			auto mu_w_s = param.water.DynamicViscosity( T_s*Xc_T, Pw_s*Xc_P, S_s );
			// neighbour
			auto mu_g_n = param.gas.DynamicViscosity( T_n*Xc_T, Pg_n*Xc_P );
			auto mu_w_n = param.water.DynamicViscosity( T_n*Xc_T, Pw_n*Xc_P, S_n );

			// 7. density
			// self
			auto rho_mass_g_s = param.gas.Density( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
			auto rho_mass_w_s = param.water.Density( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto rho_mass_h_s = param.hydrate.Density();
			// neighbour
			auto rho_mass_g_n = param.gas.Density( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
			auto rho_mass_w_n = param.water.Density( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto rho_mass_h_n = param.hydrate.Density();

			// 8. specific heat capacities
			// self
			auto Cp_g_s = param.gas.Cp( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
			auto Cp_w_s = param.water.Cp( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto rhoCp_g_s = rho_mass_g_s * Cp_g_s;
			auto rhoCp_w_s = rho_mass_w_s * Cp_w_s;
			// neighbour
			auto Cp_g_n = param.gas.Cp( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
			auto Cp_w_n = param.water.Cp( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto rhoCp_g_n = rho_mass_g_n * Cp_g_n;
			auto rhoCp_w_n = rho_mass_w_n * Cp_w_n;

			// 9. thermal conductivity
			// self
			auto kth_g_s = param.gas.ThermalConductivity( T_s*Xc_T, Pg_s*Xc_P );
			auto kth_w_s = param.water.ThermalConductivity( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto kth_h_s = param.hydrate.ThermalConductivity( T_s*Xc_T, Peff_s*Xc_P );
			auto kth_s_s = param.soil.ThermalConductivity();
			auto kth_eff_s	= (1.-porosity_s)*kth_s_s
							+ porosity_s*Sg_s*kth_g_s
							+ porosity_s*Sw_s*kth_w_s
							+ porosity_s*Sh_s*kth_h_s;
			// neighbour
			auto kth_g_n = param.gas.ThermalConductivity( T_n*Xc_T, Pg_n*Xc_P );
			auto kth_w_n = param.water.ThermalConductivity( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto kth_h_n = param.hydrate.ThermalConductivity( T_n*Xc_T, Peff_n*Xc_P );
			auto kth_s_n = param.soil.ThermalConductivity();
			auto kth_eff_n	= (1.-porosity_n)*kth_s_n
							+ porosity_n*Sg_n*kth_g_n
							+ porosity_n*Sw_n*kth_w_n
							+ porosity_n*Sh_n*kth_h_n;

			// 10. Diffusion coefficient
			// self
			auto DH2O_g_s	= tau_s * porosity_s * Sg_s * param.mixture.DiffCoeffH2OInGas( T_s*Xc_T,Pg_s*Xc_P );
			auto DCH4_w_s	= tau_s * porosity_s * Sw_s * param.mixture.DiffCoeffCH4InLiquid( T_s*Xc_T,Pw_s*Xc_P );
			auto Dc_w_s		= tau_s * porosity_s * Sw_s * param.salt.DiffCoeff( T_s*Xc_T,Pw_s*Xc_P );
			// neighbour
			auto DH2O_g_n	= tau_n * porosity_n * Sg_n * param.mixture.DiffCoeffH2OInGas( T_n*Xc_T,Pg_n*Xc_P );
			auto DCH4_w_n	= tau_n * porosity_n * Sw_n * param.mixture.DiffCoeffCH4InLiquid( T_n*Xc_T,Pw_n*Xc_P );
			auto Dc_w_n		= tau_n * porosity_n * Sw_n * param.salt.DiffCoeff( T_n*Xc_T,Pw_n*Xc_P );
			// interface values
			auto DH2O_g_int = 2.*DH2O_g_s*DH2O_g_n/( DH2O_g_s+DH2O_g_n + 1.e-15 );
			auto DCH4_w_int = 2.*DCH4_w_s*DCH4_w_n/( DCH4_w_s+DCH4_w_n + 1.e-15 );
			auto Dc_w_int 	= 2.*Dc_w_s*Dc_w_n/( Dc_w_s+Dc_w_n + 1.e-15 );


//			//upwinding wrt sediment velocity
			auto normalvelocity_s = param.parameter.SedimentVelocity()
								  * ig.unitOuterNormal(face_center_local) ;
			normalvelocity_s *= Xc_vs;

			// upwinding wrt gas-phase velocity
			auto gravity = param.parameter.g() * ig.unitOuterNormal(face_center_local) ;
			gravity *= Xc_grav;
			auto normalpotential_g = (Pg_n - Pg_s)/distance
								   + ( 0.5*(rho_mass_g_s + rho_mass_g_n) ) * gravity ;
			double omegaup_g_s = 0., omegaup_g_n = 0.;
			if( normalpotential_g>0.){
				omegaup_g_s = 0.;
				omegaup_g_n = 1.;
			}else{
				omegaup_g_s = 1.;
				omegaup_g_n = 0.;
			}
			auto normalvelocity_g = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									* ( omegaup_g_s * kr_g_s/mu_g_s + omegaup_g_n * kr_g_n/mu_g_n )
									* normalpotential_g;
			normalvelocity_g += porosity_s*Sg_s * normalvelocity_s;

			//upwinding wrt water-phase velocity
			auto normalpotential_w = (Pw_n - Pw_s)/distance
								   + ( 0.5*(rho_mass_w_s + rho_mass_w_n) ) * gravity ;
			double omegaup_w_s = 0., omegaup_w_n = 0.;
			if( normalpotential_w>0.){
				omegaup_w_s = 0.;
				omegaup_w_n = 1.;
			}else{
				omegaup_w_s = 1.;
				omegaup_w_n = 0.;
			}
			auto normalvelocity_w = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
									* ( omegaup_w_s * kr_w_s/mu_w_s + omegaup_w_n * kr_w_n/mu_w_n )
									* normalpotential_w;
			normalvelocity_w += porosity_s*Sw_s * normalvelocity_s;


			//convective flux terms
			auto convectiveflux_CH4 = ( omegaup_g_s*rho_mass_g_s*(1.-YH2O_s) + omegaup_g_n*rho_mass_g_n*(1.-YH2O_n) ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*XCH4_s + omegaup_w_n*rho_mass_w_n*XCH4_n ) * normalvelocity_w;

			auto convectiveflux_H2O = ( omegaup_g_s*rho_mass_g_s*YH2O_s + omegaup_g_n*rho_mass_g_n*YH2O_n ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*(1.-Xc_s-XCH4_s) + omegaup_w_n*rho_mass_w_n*(1.-Xc_n-XCH4_n) ) * normalvelocity_w;

			auto convectiveflux_SALT = ( omegaup_w_s * rho_mass_w_s*Xc_s + omegaup_w_n * rho_mass_w_n*Xc_n ) * normalvelocity_w;

			auto convectiveflux_HYD = porosity_s*Sh_s*rho_mass_h_s * normalvelocity_s;

			auto T_ref = param.parameter.ReferenceTemperature();
			auto convectiveflux_HEAT = (  omegaup_g_s * rhoCp_g_s * (T_s-T_ref) + omegaup_g_n * rhoCp_g_n * (T_n-T_ref)  ) * normalvelocity_g
									 + (  omegaup_w_s * rhoCp_w_s * (T_s-T_ref) + omegaup_w_n * rhoCp_w_n * (T_n-T_ref)  ) * normalvelocity_w;

			//diffusive flux terms
			double j_CH4_w  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DCH4_w_int * ( XCH4_n - XCH4_s )/distance;
			double j_SALT_w = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * Dc_w_int   * ( Xc_n   - Xc_s   )/distance;
			double j_H2O_w  = - ( j_CH4_w + j_SALT_w );
			double j_H2O_g  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DH2O_g_int * ( YH2O_n - YH2O_s )/distance;;
			double j_CH4_g  = - j_H2O_g;

			auto diffusiveflux_CH4 = j_CH4_w + j_CH4_g;

			auto diffusiveflux_H2O = j_H2O_w + j_H2O_g;

			auto diffusiveflux_SALT = j_SALT_w;

			auto diffusiveflux_HYD = 0.;

			auto diffusiveflux_HEAT = - ( 2.*kth_eff_s*kth_eff_n/( kth_eff_s+kth_eff_n ) ) * ( T_n - T_s )/distance;


			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4 ;
			r_s.accumulate(lfs_Pw_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Pw_n, 0, -tmp*face_volume);

			// H2O-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O ;
			r_s.accumulate(lfs_Sg_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Sg_n, 0, -tmp*face_volume);

			// SALT-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;
			r_s.accumulate(lfs_Xc_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Xc_n, 0, -tmp*face_volume);

			// HYD-phase-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_HYD + Xc_diff_m * diffusiveflux_HYD ;
			r_s.accumulate(lfs_Sh_s, 0, +tmp*face_volume);
			r_n.accumulate(lfs_Sh_n, 0, -tmp*face_volume);

			// ENERGY balance
			tmp = Xc_conv_h * convectiveflux_HEAT + Xc_diff_h * diffusiveflux_HEAT ;
			r_s.accumulate(lfs_T_s , 0, +tmp*face_volume);
			r_n.accumulate(lfs_T_n , 0, -tmp*face_volume);

	  }

	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_boundary ( const IG& ig,
			  	  	  	  	const LFSU& lfsu, const X& x, const LFSV& lfsv,
							R& r ) const
	  {
			typedef typename LFSU::template Child<Indices::PVId_Pw>::Type LFS_Pw_s;
			const LFS_Pw_s& lfs_Pw_s = lfsu.template child<Indices::PVId_Pw>();
			typedef typename LFSU::template Child<Indices::PVId_Sg>::Type LFS_Sg_s;
			const LFS_Sg_s& lfs_Sg_s = lfsu.template child<Indices::PVId_Sg>();
			typedef typename LFSU::template Child<Indices::PVId_XCH4>::Type LFS_XCH4_s;
			const LFS_XCH4_s& lfs_XCH4_s = lfsu.template child<Indices::PVId_XCH4>();
			typedef typename LFSU::template Child<Indices::PVId_YH2O>::Type LFS_YH2O_s;
			const LFS_YH2O_s& lfs_YH2O_s = lfsu.template child<Indices::PVId_YH2O>();
			typedef typename LFSU::template Child<Indices::PVId_Xc>::Type LFS_Xc_s;
			const LFS_Xc_s& lfs_Xc_s = lfsu.template child<Indices::PVId_Xc>();
			typedef typename LFSU::template Child<Indices::PVId_Sh>::Type LFS_Sh_s;
			const LFS_Sh_s& lfs_Sh_s = lfsu.template child<Indices::PVId_Sh>();
			typedef typename LFSU::template Child<Indices::PVId_T >::Type LFS_T_s ;
			const LFS_T_s&  lfs_T_s  = lfsu.template child<Indices::PVId_T >();

	        // References to inside and outside cells
	        const auto& cell_inside  = ig.inside();
			const IndexSet &indexSet = gv.indexSet();
			int inside_cell_number  = indexSet.index(cell_inside);

	        // get geometries
	        auto geo = ig.geometry();
	        const auto dim = geo.mydimension;
	        auto geo_inside = cell_inside.geometry();

	        // cell geometries
	        auto ref_el_inside 	= referenceElement(geo_inside);
	        auto inside_cell_center_local 	= ref_el_inside.position(0,0);
	        auto inside_cell_center_global 	= geo_inside.center();

	        // face geometry
	        auto ref_el = referenceElement(geo);
	        auto face_center_local = ref_el.position(0,0);
	        auto face_center_global = geo.center();
	        auto face_volume = geo.volume();
	        auto normal = ig.unitOuterNormal(face_center_local);
	        face_volume *= param.mesh.areaFactor( inside_cell_center_global[0],
	        									  face_center_global[0],
												  normal);

	        // distance of cell centers
	        auto d = geo.global(face_center_local);
	        d -= inside_cell_center_global;
	        auto distance = d.two_norm();

			// compute primary vars at local self and neighbour centers
	        double Pw_s   = x(lfs_Pw_s,  0);
			double Sg_s   = x(lfs_Sg_s,  0);
			double XCH4_s = x(lfs_XCH4_s,0);
			double YH2O_s = x(lfs_YH2O_s,0);
			double Xc_s   = x(lfs_Xc_s,  0);
			double Sh_s   = x(lfs_Sh_s,  0);
			double T_s    = x(lfs_T_s ,  0);

			double Pw_n   = Pw_s;
			double Sg_n   = Sg_s;
			double XCH4_n = XCH4_s;
			double YH2O_n = YH2O_s;
			double Xc_n   = Xc_s;
			double Sh_n   = Sh_s;
			double T_n    = T_s;

			// evaluate boundary condition types for {Pw,Sg} or {Fw,Fg}
			auto bctype = bc.type(ig,face_center_local) ;

			// evaluate boundary condition values for {Pw,Sg} or {Fw,Fg}
			auto bcvalue = bc.value(ig,face_center_local) ;

			/* SET DIRICHLET VALUES FOR FACE (PSEUDO-NEIGHBOUR) FOR Pw, Sg, XCH4, Xc, and T
			 * note:- For setting dirichlet value for XCH4, we need to calculate the equilibrium conc. first...
			 * bcvalue_XCH4 and bcvalue_YH2O are assigned after 4.secondary variables */
			if( bctype[Indices::BCId_water] == Indices::dirichlet ){
				Pw_n = bcvalue[Indices::BCId_water];
			}

			if( bctype[Indices::BCId_gas] == Indices::dirichlet ){
				Sg_n = bcvalue[Indices::BCId_gas];
			}

			if( bctype[Indices::BCId_salt] == Indices::dirichlet ){
				Xc_n = bcvalue[Indices::BCId_salt];
			}

			if( bctype[Indices::BCId_heat]  == Indices::dirichlet ){
				T_n  = bcvalue[Indices::BCId_heat] ;
			}


			// compute properties at local cell center and neighbour cell center
			// Hydraulic properties for self and neighbour
			// 1. porosity
			auto porosity_s = param.soil.SedimentPorosity( cell_inside,inside_cell_center_local );
#ifdef COMPACTION
			porosity_s *= param.parameter.CompactionFunction( inside_cell_center_global );
#endif
			auto porosity_n = porosity_s;
			// 2. abs. permeability
			auto permeability_s =  param.soil.SedimentPermeabilityVector( cell_inside,inside_cell_center_local ) * ig.unitOuterNormal(face_center_local);
			permeability_s *= param.hydraulicProperty.PermeabilityScalingFactor( cell_inside,inside_cell_center_local, Sh_s, porosity_s );
			permeability_s = std::abs(permeability_s);
			auto permeability_n = permeability_s;
			// 3. tortuosity
			auto tau_s = param.soil.Tortuosity( porosity_s*(1.-Sh_s) );
			auto tau_n = param.soil.Tortuosity( porosity_n*(1.-Sh_n) );

			// 4. secondary variables
			// at self
			auto S_s = Xc_s * (param.salt.MolarMass()/param.water.MolarMass());
			auto Sw_s = 1.-Sh_s-Sg_s;
			auto Pc_s = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_s, Sh_s, porosity_s );
			auto Pg_s = Pw_s + Pc_s;
			auto Peff_s = ( Pg_s*Sg_s + Pw_s*Sw_s )/(1.-Sh_s);
			auto zCH4_s = param.eos.EvaluateCompressibilityFactor( T_s*Xc_T,Pg_s*Xc_P );
			// at neighbour
			auto S_n = Xc_n * (param.salt.MolarMass()/param.water.MolarMass());
			auto Sw_n = 1.-Sh_n-Sg_n;
			auto Pc_n = param.hydraulicProperty.CapillaryPressure( cell_inside,inside_cell_center_local, Sw_n, Sh_n, porosity_n );
			auto Pg_n = Pw_n + Pc_n;
			auto Peff_n = ( Pg_n*Sg_n + Pw_n*Sw_n )/(1.-Sh_n);
			auto zCH4_n = param.eos.EvaluateCompressibilityFactor( T_n*Xc_T,Pg_n*Xc_P );

			// 5. relative permeability
			// self
			auto kr_g_s = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_s, Sh_s );
			auto kr_w_s = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_s, Sh_s );
			// neighbour
			auto kr_g_n = param.hydraulicProperty.krg( cell_inside, inside_cell_center_local, Sw_n, Sh_n );
			auto kr_w_n = param.hydraulicProperty.krw( cell_inside, inside_cell_center_local, Sw_n, Sh_n );

			// 6. dynamic viscosity
			// self
			auto mu_g_s = param.gas.DynamicViscosity( T_s*Xc_T, Pg_s*Xc_P );
			auto mu_w_s = param.water.DynamicViscosity( T_s*Xc_T, Pw_s*Xc_P, S_s );
			// neighbour
			auto mu_g_n = param.gas.DynamicViscosity( T_n*Xc_T, Pg_n*Xc_P );
			auto mu_w_n = param.water.DynamicViscosity( T_n*Xc_T, Pw_n*Xc_P, S_n );

			// 7. density
			// self
			auto rho_mass_g_s = param.gas.Density( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
			auto rho_mass_w_s = param.water.Density( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto rho_mass_h_s = param.hydrate.Density();
			// neighbour
			auto rho_mass_g_n = param.gas.Density( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
			auto rho_mass_w_n = param.water.Density( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto rho_mass_h_n = param.hydrate.Density();

			// 8. specific heat capacities
			// self
			auto Cp_g_s = param.gas.Cp( T_s*Xc_T, Pg_s*Xc_P, zCH4_s );
			auto Cp_w_s = param.water.Cp( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto rhoCp_g_s = rho_mass_g_s * Cp_g_s;
			auto rhoCp_w_s = rho_mass_w_s * Cp_w_s;
			// neighbour
			auto Cp_g_n = param.gas.Cp( T_n*Xc_T, Pg_n*Xc_P, zCH4_n );
			auto Cp_w_n = param.water.Cp( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto rhoCp_g_n = rho_mass_g_n * Cp_g_n;
			auto rhoCp_w_n = rho_mass_w_n * Cp_w_n;

			// 9. thermal conductivity
			// self
			auto kth_g_s = param.gas.ThermalConductivity( T_s*Xc_T, Pg_s*Xc_P );
			auto kth_w_s = param.water.ThermalConductivity( T_s*Xc_T, Pw_s*Xc_P, S_s );
			auto kth_h_s = param.hydrate.ThermalConductivity( T_s*Xc_T, Peff_s*Xc_P );
			auto kth_s_s = param.soil.ThermalConductivity();
			auto kth_eff_s	= (1.-porosity_s)*kth_s_s
							+ porosity_s*Sg_s*kth_g_s
							+ porosity_s*Sw_s*kth_w_s
							+ porosity_s*Sh_s*kth_h_s;
			// neighbour
			auto kth_g_n = param.gas.ThermalConductivity( T_n*Xc_T, Pg_n*Xc_P );
			auto kth_w_n = param.water.ThermalConductivity( T_n*Xc_T, Pw_n*Xc_P, S_n );
			auto kth_h_n = param.hydrate.ThermalConductivity( T_n*Xc_T, Peff_n*Xc_P );
			auto kth_s_n = param.soil.ThermalConductivity();
			auto kth_eff_n	= (1.-porosity_n)*kth_s_n
							+ porosity_n*Sg_n*kth_g_n
							+ porosity_n*Sw_n*kth_w_n
							+ porosity_n*Sh_n*kth_h_n;

			// 10. Diffusion coefficient
			// self
			auto DH2O_g_s	= tau_s * porosity_s * Sg_s * param.mixture.DiffCoeffH2OInGas( T_s*Xc_T,Pg_s*Xc_P );
			auto DCH4_w_s	= tau_s * porosity_s * Sw_s * param.mixture.DiffCoeffCH4InLiquid( T_s*Xc_T,Pw_s*Xc_P );
			auto Dc_w_s		= tau_s * porosity_s * Sw_s * param.salt.DiffCoeff( T_s*Xc_T,Pw_s*Xc_P );
			// neighbour
			auto DH2O_g_n	= tau_n * porosity_n * Sg_n * param.mixture.DiffCoeffH2OInGas( T_n*Xc_T,Pg_n*Xc_P );
			auto DCH4_w_n	= tau_n * porosity_n * Sw_n * param.mixture.DiffCoeffCH4InLiquid( T_n*Xc_T,Pw_n*Xc_P );
			auto Dc_w_n		= tau_n * porosity_n * Sw_n * param.salt.DiffCoeff( T_n*Xc_T,Pw_n*Xc_P );
			// interface values
			auto DH2O_g_int = 2.*DH2O_g_s*DH2O_g_n/( DH2O_g_s+DH2O_g_n + 1.e-15 );
			auto DCH4_w_int = 2.*DCH4_w_s*DCH4_w_n/( DCH4_w_s+DCH4_w_n + 1.e-15 );
			auto Dc_w_int 	= 2.*Dc_w_s*Dc_w_n/( Dc_w_s+Dc_w_n + 1.e-15 );


//			//upwinding wrt sediment velocity
			auto normalvelocity_s = 0.;//
			normalvelocity_s *= param.parameter.SedimentVelocity()
								  * ig.unitOuterNormal(face_center_local) ;
			normalvelocity_s *= Xc_vs ;
//			double omegaup_s_s = 0., omegaup_s_n = 0.;
//			if( normalvelocity_s>0.){
//				omegaup_s_s = 0.;
//				omegaup_s_n = 1.;
//			}else{
//				omegaup_s_s = 1.;
//				omegaup_s_n = 0.;
//			}

			// upwinding wrt gas-phase velocity
			auto gravity = param.parameter.g() * ig.unitOuterNormal(face_center_local) ;
			gravity *= Xc_grav;
			auto normalpotential_g = (Pg_n - Pg_s)/distance
								   + ( 0.5*(rho_mass_g_s + rho_mass_g_n) ) * gravity ;
			double omegaup_g_s = 0., omegaup_g_n = 0.;
			if( normalpotential_g>0.){
				omegaup_g_s = 0.;
				omegaup_g_n = 1.;
			}else{
				omegaup_g_s = 1.;
				omegaup_g_n = 0.;
			}

			//upwinding wrt water-phase velocity
			auto normalpotential_w = (Pw_n - Pw_s)/distance
								   + ( 0.5*(rho_mass_w_s + rho_mass_w_n) ) * gravity ;
			double omegaup_w_s = 0., omegaup_w_n = 0.;
			if( normalpotential_w>0.){
				omegaup_w_s = 0.;
				omegaup_w_n = 1.;
			}else{
				omegaup_w_s = 1.;
				omegaup_w_n = 0.;
			}

			double normalvelocity_g =0., normalvelocity_w =0.;
			if( bctype[Indices::BCId_water] == Indices::neumann and ( bctype[Indices::BCId_gas] == Indices::neumann ) ){
				normalvelocity_g = bcvalue[Indices::BCId_gas];
				normalvelocity_w = bcvalue[Indices::BCId_water];
			}
			else if( bctype[Indices::BCId_water] == Indices::dirichlet and ( bctype[Indices::BCId_gas] == Indices::neumann ) ){
				normalvelocity_g = bcvalue[Indices::BCId_gas];
				normalvelocity_w = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
								   * ( omegaup_w_s * kr_w_s/mu_w_s + omegaup_w_n * kr_w_n/mu_w_n )
								   * normalpotential_w;
			}
			else if( bctype[Indices::BCId_water] == Indices::neumann and ( bctype[Indices::BCId_gas] == Indices::dirichlet ) ){
				normalvelocity_g = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
								   * ( omegaup_g_s * kr_g_s/mu_g_s + omegaup_g_n * kr_g_n/mu_g_n )
								   * normalpotential_g;
				normalvelocity_w = bcvalue[Indices::BCId_water];
			}
			else if( bctype[Indices::BCId_water] == Indices::dirichlet and ( bctype[Indices::BCId_gas] == Indices::dirichlet ) ){
				normalvelocity_g = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
								   * ( omegaup_g_s * kr_g_s/mu_g_s + omegaup_g_n * kr_g_n/mu_g_n )
								   * normalpotential_g;
				normalvelocity_w = - ( 2.*permeability_s*permeability_n/(permeability_s+permeability_n) )
								   * ( omegaup_w_s * kr_w_s/mu_w_s + omegaup_w_n * kr_w_n/mu_w_n )
								   * normalpotential_w;
			}
			normalvelocity_g += porosity_s*Sg_s * normalvelocity_s;
			normalvelocity_w += porosity_s*Sw_s * normalvelocity_s;

			double temperature_gradient = 0.;
			if( bctype[Indices::BCId_heat] == Indices::neumann ){
				temperature_gradient = bcvalue[Indices::BCId_heat];
			}else{
				temperature_gradient = (T_n-T_s)/distance;
			}

			/* FLUX TERMS */

			//convective flux terms
			auto convectiveflux_CH4 = ( omegaup_g_s*rho_mass_g_s*(1.-YH2O_s) + omegaup_g_n*rho_mass_g_n*(1.-YH2O_n) ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*XCH4_s + omegaup_w_n*rho_mass_w_n*XCH4_n ) * normalvelocity_w;

			auto convectiveflux_H2O = ( omegaup_g_s*rho_mass_g_s*YH2O_s + omegaup_g_n*rho_mass_g_n*YH2O_n ) * normalvelocity_g
									+ ( omegaup_w_s*rho_mass_w_s*(1.-Xc_s-XCH4_s) + omegaup_w_n*rho_mass_w_n*(1.-Xc_n-XCH4_n) ) * normalvelocity_w;

			auto convectiveflux_SALT = ( omegaup_w_s * rho_mass_w_s*Xc_s + omegaup_w_n * rho_mass_w_n*Xc_n ) * normalvelocity_w;

			auto convectiveflux_HYD = porosity_s*Sh_s*rho_mass_h_s * normalvelocity_s;

			auto T_ref = param.parameter.ReferenceTemperature();
			auto convectiveflux_HEAT = (  omegaup_g_s * rhoCp_g_s * (T_s-T_ref) + omegaup_g_n * rhoCp_g_n * (T_n-T_ref)  ) * normalvelocity_g
									 + (  omegaup_w_s * rhoCp_w_s * (T_s-T_ref) + omegaup_w_n * rhoCp_w_n * (T_n-T_ref)  ) * normalvelocity_w;

			//diffusive flux terms
			double j_CH4_w  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DCH4_w_int * ( XCH4_n - XCH4_s )/distance;
			double j_SALT_w = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * Dc_w_int   * ( Xc_n   - Xc_s   )/distance;
			double j_H2O_w  = - ( j_CH4_w + j_SALT_w );
			double j_H2O_g  = - 0.5*( rho_mass_w_s + rho_mass_w_n ) * DH2O_g_int * ( YH2O_n - YH2O_s )/distance;;
			double j_CH4_g  = - j_H2O_g;

			auto diffusiveflux_CH4 = j_CH4_w + j_CH4_g;

			auto diffusiveflux_H2O = j_H2O_w + j_H2O_g;

			auto diffusiveflux_SALT = j_SALT_w;

			auto diffusiveflux_HYD = 0.;

			auto diffusiveflux_HEAT = - ( 2.*kth_eff_s*kth_eff_n/( kth_eff_s+kth_eff_n ) ) * temperature_gradient;


			/*ACCCUMULATE RESIDUALS*/
			double tmp=0.;

			// CH4-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_CH4 + Xc_diff_m * diffusiveflux_CH4 ;
			r.accumulate(lfs_Pw_s, 0, +tmp*face_volume);

			// H2O-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_H2O + Xc_diff_m * diffusiveflux_H2O ;
			r.accumulate(lfs_Sg_s, 0, +tmp*face_volume);

			// SALT-component-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_SALT + Xc_diff_m * diffusiveflux_SALT ;
			r.accumulate(lfs_Xc_s, 0, +tmp*face_volume);

			// HYD-phase-wise mass-balance
			tmp = Xc_conv_m * convectiveflux_HYD + Xc_diff_m * diffusiveflux_HYD ;
			r.accumulate(lfs_Sh_s, 0, +tmp*face_volume);

			// ENERGY balance
			tmp = Xc_conv_h * convectiveflux_HEAT + Xc_diff_h * diffusiveflux_HEAT ;
			r.accumulate(lfs_T_s , 0, +tmp*face_volume);
	  }

};
