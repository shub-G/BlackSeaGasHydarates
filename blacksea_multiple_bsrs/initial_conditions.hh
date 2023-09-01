template<typename GV,typename Properties,typename PGMap>
class ProblemInitialConditions
{
private:
	  const GV& gv;
	  const Properties& property;
	  const PGMap& pg ;
	  const static int dim = GV::dimension;
	  constexpr static double eps = 1.e-6;

public:

	  //! construct from grid view
	  ProblemInitialConditions ( const GV& gv_,
			  	  	  	  	 	 const Properties& property_,
								 const PGMap& pg_ )
	  : gv( gv_ ) ,
		property(property_),
		pg( pg_ )
	  {}

	  /* Initial Conditions */
	  std::vector< double >
	  evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
			  	const Dune::FieldVector<double,dim>& xlocal) const {

		  std::vector< double > icvalue(Indices::numOfPVs,0.);

		  Dune::FieldVector<double,dim> xglobal/*ndim*/ = element.geometry().global(xlocal);

		  double x_c = property.characteristicValue.x_c;
		  double Z_max /*ndim*/ = property.mesh.z_PSFC(xglobal[0]/*ndim*/);
		  /******************************************************************************/
		  // SATURATIONS
		  double Sg = 0.0;
		  double Sh = property.parameter.InitialHydrateSaturation(xglobal);
		  double Sw = 1.-Sg-Sh;
		  icvalue[Indices::PVId_Sg] = Sg ;
		  icvalue[Indices::PVId_Sh] = Sh ;

		  /******************************************************************************/
		  // PRESSURES
		  double porosity = property.soil.SedimentPorosity(element,xlocal);
		  double Pc = property.hydraulicProperty.CapillaryPressure(element,xlocal,Sw,Sh,porosity)
				    * property.characteristicValue.P_c; /*Pa*/
		  double Pw_top = property.parameter.PSF_Pressure(); /*Pa*/
		  double Pw = Pw_top + 1000.*10.*(Z_max-xglobal[1])*x_c; /*Pa*/
		  double Pg = Pw + Pc; /*Pa*/
		  icvalue[Indices::PVId_Pw] = Pw/property.characteristicValue.P_c ;/*ndim*/

		  /******************************************************************************/
		  // TEMPERATURE
		  double T_top 	= property.parameter.PSF_Temperature(); /*K*/
		  double T_grad	= property.parameter.RegionalThermalGradient(); /*degC per m*/
		  double T = T_top;
		  if(xglobal[1]<(Z_max)+eps){
			  T +=  T_grad*(Z_max-xglobal[1])*x_c ;
		  }
		  icvalue[Indices::PVId_T] = T/property.characteristicValue.T_c;

		  /******************************************************************************/
		  // SALINITY
		  double sal = property.parameter.PaleoSalinity();
		  double xc = sal * (property.water.MolarMass()/property.salt.MolarMass()) ;
		  icvalue[Indices::PVId_Xc ] = xc;

		  /******************************************************************************/
		  // MOLE FRACTIONS
		  auto zCH4 = property.eos.EvaluateCompressibilityFactor( T,Pg );
		  auto Xf = property.mixture.EquilibriumMoleFractions( T,Pg, xc, zCH4);

		  auto XCH4_eqb = Xf[Indices::compId_XCH4];
		  auto nxch4 = property.parameter.InitialXCH4Fraction();
		  double xch4 = nxch4*XCH4_eqb;
		  icvalue[Indices::PVId_XCH4 ] = xch4;

		  auto YH2O_eqb = Xf[Indices::compId_YH2O];
		  double yh2o = YH2O_eqb;
		  icvalue[Indices::PVId_YH2O ] = yh2o;

		  /******************************************************************************/

		  return icvalue; /*ndim*/
	  }

	  //! get a reference to the grid view
	  inline const GV& getGridView () {return gv;}
};
