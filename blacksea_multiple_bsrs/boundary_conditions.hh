template<typename GV, typename Properties, typename PGMap >
class ProblemBoundaryConditions
{
private :
	const GV& gv ;
	const Properties& property;
	const PGMap& pgmap ;

	const static int dim = GV::dimension;

	ProblemInitialConditions<GV,Properties,PGMap> icvalue;
	constexpr static double eps = 1.e-3;

public:

	// ! construct from gridview
	ProblemBoundaryConditions ( const GV& gv_, 
								const Properties& property_, 
								const PGMap& pgmap_ )
	: gv ( gv_ ) ,
	  property(property_),
	  pgmap ( pgmap_ ),
	  icvalue(gv_,property_,pgmap_)
	{}

	/* boundary types */
	template<typename I> std::vector< int >
	type( I& intersection,
		  const Dune::FieldVector<double,dim-1>& xlocal ) const {

		auto xglobal = intersection.geometry().global(xlocal);

		std::vector< int > bct(Indices::numOfBCs,Indices::neumann);
		if( property.mesh.isTopBoundary(xglobal) ){
			bct[Indices::BCId_water] = Indices::dirichlet;
			bct[Indices::BCId_gas] 	 = Indices::dirichlet;
			bct[Indices::BCId_salt ] = Indices::dirichlet;
			bct[Indices::BCId_heat ] = Indices::dirichlet;
		}

		return bct;
	}

	/* boundary values */
	template<typename I> std::vector<double>
	value ( I& intersection,
		    const Dune::FieldVector<double,dim-1>& xlocal ) const
	{
		std::vector< double > bcvalue(Indices::numOfBCs,0.); // DEFAULT is NEUMANN 0

		auto xglobal = intersection.geometry().global(xlocal);
		const auto& element = intersection.inside();
        auto element_xlocal = referenceElement(element.geometry()).position(0,0);

        double sedimentation_depth /*m*/  = property.parameter.SedimentationDepth( xglobal[0]*property.characteristicValue.x_c );
        double grad_T = property.parameter.RegionalThermalGradient();//deg C/m

		if( property.mesh.isTopBoundary(xglobal) ){
			auto icv /*ndim*/ = icvalue.evaluate(element,element_xlocal);

			double porosity_top = property.soil.SedimentPorosity(element,element_xlocal);
			double Pw_top = icv[Indices::PVId_Pw]
							+ 1000.*property.parameter.g()[dim-1]*sedimentation_depth/property.characteristicValue.P_c ;
			double Sg_top = icv[Indices::PVId_Sg];
			double xc_top = icv[Indices::PVId_Xc];
			double T_top  = icv[Indices::PVId_T]
							+ grad_T * sedimentation_depth /property.characteristicValue.T_c;

			bcvalue[Indices::BCId_water] = Pw_top ;
			bcvalue[Indices::BCId_gas  ] = Sg_top ;
			bcvalue[Indices::BCId_salt ] = xc_top ;
			bcvalue[Indices::BCId_heat ] = T_top  ;

		}
		else if( property.mesh.isBottomBoundary(xglobal) ){
			bcvalue[Indices::BCId_heat ] = grad_T* (property.characteristicValue.x_c/property.characteristicValue.T_c);
		}

		return bcvalue /*ndim*/;
	}

};



