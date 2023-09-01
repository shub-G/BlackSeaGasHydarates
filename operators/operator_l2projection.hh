template <class GV, class PARAMS, class GFS, class U>
class LocalOperatorPROJ :
	public Dune::PDELab::NumericalJacobianApplyVolume		<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::NumericalJacobianVolume			<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::NumericalJacobianApplySkeleton		<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::NumericalJacobianSkeleton			<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::NumericalJacobianApplyBoundary		<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::NumericalJacobianBoundary			<LocalOperatorPROJ<GV,PARAMS,GFS,U> >,
	public Dune::PDELab::FullSkeletonPattern,                     // matrix entries skeleton
	public Dune::PDELab::FullVolumePattern,
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
		    const GV& gv;
		    const PARAMS& param;
			GFS gfs;
			U *ufv;
			unsigned int intorder ;
			double epsilon;
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
			enum { doPatternVolume 	 = true  };
			enum { doPatternSkeleton = false };

			// residual assembly flags
			enum { doAlphaVolume  	= true 	};
			enum { doLambdaVolume  	= true 	};
			enum { doAlphaSkeleton  = false };
			enum { doAlphaBoundary  = false };
			enum { doLambdaBoundary = false };

			typedef typename GV::IndexSet IndexSet;

			typedef Dune::PDELab::LocalFunctionSpace<GFS> LFS;
		    typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
			typedef typename U::template LocalView<LFSCache> VectorView;

			// constructor stores parameters
			LocalOperatorPROJ(	const GV& gv_,
								const PARAMS& param_,
								GFS 	gfs_,
								U		*ufv_,
								unsigned int 	intorder_ = 4,
								double 			epsilon_ = 1.e-6)
			: gv(gv_),
			  param(param_),
			  gfs(gfs_),
			  ufv(ufv_),
			  intorder( intorder_ ),
			  epsilon( epsilon_ )
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

			// volume integral depending on test and ansatz functions
			template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
			void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
			{
				// define types
				using RF = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
				using RangeType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
				using JacobianType = typename LFSU::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

				// dimensions
				const int dim = EG::Entity::dimension;

				// Get cell
				const auto& cell = eg.entity();

				// Get geometry
				auto geo = eg.geometry();

				// evaluate diffusion tensor at cell center, assume it is constant over elements
				auto ref_el = referenceElement(geo);
				auto localcenter = ref_el.position(0,0);

				// Transformation matrix
				typename EG::Geometry::JacobianInverseTransposed jac;

				// loop over quadrature points
				for (const auto& ip : quadratureRule(geo,intorder))
				{
					// evaluate basis functions
					std::vector< std::vector<RangeType> > phi( Indices::numOfPhases );
			        for( int i=0; i<phi.size(); i++ ){
			        	 phi[i] = std::vector<RangeType> ( lfsu.child(i).size() );
			        	 lfsu.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),phi[i]);
			        }

					std::vector< std::vector<RangeType> > psi( Indices::numOfPhases );
			        for( int i=0; i<psi.size(); i++ ){
			        	 psi[i] = std::vector<RangeType> ( lfsv.child(i).size() );
			        	 lfsv.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),psi[i]);
			        }

					auto ip_global = geo.global(ip.position());

					// evaluate u (u->phase pressures)
			        std::vector<double> u( Indices::numOfPhases ,0.);
			        for(int i=0; i<u.size(); i++ ){
			        	for (int j=0; j<lfsu.child(i).size(); j++)
			        		u[i] += x(lfsu.child(i),j) * phi[i][j];
			        }

					// evaluate gradient of basis functions
			        std::vector< std::vector<JacobianType> > jsu( Indices::numOfPhases );
			        for( int i=0; i<jsu.size(); i++ ){
			        	 jsu[i] = std::vector<JacobianType> ( lfsu.child(i).size() );
			        	 lfsu.child(i).finiteElement().localBasis().evaluateJacobian(ip.position(),jsu[i]);
			        }

			        std::vector< std::vector<JacobianType> > jsv( Indices::numOfPhases );
			        for( int i=0; i<jsv.size(); i++ ){
			        	 jsv[i] = std::vector<JacobianType> ( lfsv.child(i).size() );
			        	 lfsv.child(i).finiteElement().localBasis().evaluateJacobian(ip.position(),jsv[i]);
			        }

					// transform gradients of shape functions to real element
					jac = geo.jacobianInverseTransposed(ip.position());

					// evaluade gradients of shape fncs.
			        std::vector< std::vector < Dune::FieldVector<RF,dim > > > gradphi( Indices::numOfPhases );
			        for( int i=0; i<gradphi.size(); i++ ){
			        	gradphi[i] = std::vector<Dune::FieldVector<RF,dim>> ( lfsu.child(i).size() );
			        	for (int j=0; j<lfsu.child(i).size(); j++){
			        	          jac.mv(jsu[i][j][0],gradphi[i][j]);
			        	}
			        }

			        std::vector< std::vector < Dune::FieldVector<RF,dim > > > gradpsi( Indices::numOfPhases );
			        for( int i=0; i<gradpsi.size(); i++ ){
			        	gradpsi[i] = std::vector<Dune::FieldVector<RF,dim>> ( lfsv.child(i).size() );
			        	for (int j=0; j<lfsv.child(i).size(); j++){
			        	          jac.mv(jsv[i][j][0],gradpsi[i][j]);
			        	}
			        }

					// compute gradient of u (phase pressures)
			        std::vector< Dune::FieldVector<RF,dim > > gradu( Indices::numOfPhases );
					for( int i=0; i<gradu.size(); i++ ){
						gradu[i] =  Dune::FieldVector<RF,dim>(0.0);
						for (int j=0; j<lfsu.child(i).size(); j++){
							gradu[i].axpy(x(lfsu.child(i),j),gradphi[i][j]);
						}
					}

					// integrate
					// FV --> FEM
					//|| u_FV - u_FE || --> Min
					//
					RF factor = ip.weight() * geo.integrationElement(ip.position());
					for (int i=0; i<Indices::numOfPhases; i++){
						double tmp = 0.;
						for (int j=0; j<lfsv.child(i).size(); j++){
							tmp = epsilon*( gradu[i]*gradpsi[i][j] ) + u[i] * psi[i][j] ;
							r.accumulate(lfsv.child(i),j,tmp*factor );
						}
					}

				}//End Quadrature Rule

			}


			// volume integral depending only on test functions
			template<typename EG, typename LFSV, typename R>
			void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
			{
				// define types
				using RF = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
				using RangeType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
				using JacobianType = typename LFSV::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;

				// dimensions
				const int dim = EG::Entity::dimension;

				// Get cell
				const auto& cell = eg.entity();

				// Get geometry
				auto geo = eg.geometry();

				// evaluate diffusion tensor at cell center, assume it is constant over elements
				auto ref_el = referenceElement(geo);
				auto localcenter = ref_el.position(0,0);
		        auto globalcenter = cell.geometry().global(localcenter);

				LFS lfs( gfs );
				LFSCache lfscache( lfs );
				VectorView ufv_view( (*ufv) );

				lfs.bind( (cell) );
				lfscache.update();
				ufv_view.bind( lfscache );

				std::vector<RF> ul_fv(lfs.size());
				ufv_view.read( ul_fv );

				// evaluate right hand side parameter function
				auto Pw_fv = ul_fv[lfs.child(Indices::PVId_Pw).localIndex(0)];
				auto Sg_fv = ul_fv[lfs.child(Indices::PVId_Sg).localIndex(0)];
				auto Sh_fv = ul_fv[lfs.child(Indices::PVId_Sh).localIndex(0)];

				auto Sw_fv = 1.-Sg_fv-Sh_fv;
				auto porosity = param.soil.SedimentPorosity( cell,localcenter );
	#ifdef COMPACTION
				porosity *= param.parameter.CompactionFunction( globalcenter );
	#endif
				auto Pc_fv = param.hydraulicProperty.CapillaryPressure( cell,localcenter, Sw_fv, Sh_fv, porosity );
				auto Pg_fv = Pw_fv + Pc_fv;

				// loop over quadrature points
				for (const auto& ip : quadratureRule(geo,intorder))
				{
					// evaluate shape functions
					std::vector<std::vector<RangeType>> psi( Indices::numOfPhases );
			        for( int i=0; i<psi.size(); i++ ){
			        	 psi[i] = std::vector<RangeType> ( lfsv.child(i).size() );
			        	 lfsv.child(i).finiteElement().localBasis().evaluateFunction(ip.position(),psi[i]);
			        }

					// integrate
					RF factor = ip.weight() * geo.integrationElement(ip.position());

					for (int j=0; j<lfsv.child(Indices::phaseId_g).size(); j++){
						r.accumulate(lfsv.child(Indices::phaseId_g),j,( -Pg_fv * psi[Indices::phaseId_g][j] )*factor);
					}

					for (int j=0; j<lfsv.child(Indices::phaseId_w).size(); j++){
						r.accumulate(lfsv.child(Indices::phaseId_w),j,( -Pw_fv * psi[Indices::phaseId_w][j] )*factor);
					}

				}//END: quadrature rule

			}//END: lambda_volume

	  };

