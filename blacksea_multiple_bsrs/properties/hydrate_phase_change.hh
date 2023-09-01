template<typename GV, typename Parameters, typename PGMap>
class HydratePhaseChangeKinetics
{
private:
	  const GV& gv;
	  const Parameters& parameter;
	  const PGMap& pgmap ;

	  const static int dim = GV::dimension;

	  CharacteristicValues characteristicValue;
	  Methane methane;
	  Water water;
	  Hydrate hydrate;
	  Salt salt;
	  Soil<GV,Parameters,PGMap > soil;
	  HydraulicProperties<GV,Parameters,PGMap> hydraulicProperty;

public:

  //! construct from grid view
  HydratePhaseChangeKinetics ( const GV& gv_ , 
		  	  	  	  	  	   const Parameters& parameter_ , 
							   const PGMap& pgmap_ )
  : gv( gv_ ), 
	parameter(parameter_), 
	pgmap( pgmap_ ),
	soil(gv_,parameter_,pgmap_),
	hydraulicProperty(gv_,parameter_,pgmap_)
  {}

	// EQULIBRIUM CONDITIONS FOR HYDRATE:
	//------------------------------------
	//equilibrium pressure
	double EquilibriumPressure (double T/*K*/,double S) const {

//		double A = 38.592 , 	B = 8533.8 	,	C = 16.32;
//		double P_eq = 1.e3 * exp( A - B/( T ) + C*S ); // defined in Pascals

//		double P_eq = 1.0e6*exp(- 1.6892692e3 - 0.15162984*T + 5.60912482e4/T + 2.72067506E2*log(T)
//								+ (S/1.80655e-3)*( 2.06621298e4 + 14.18127*T - 8.3417066e-3*T*T - 3.78757519E5/T - 4.01544417e3*log(T) )
//								+ (S/1.80655e-3)*(S/1.80655e-3)*( 6.04018045e2 + 0.420253434*T - 2.501548e-4*T*T - 1.09745335E4/T - 1.17640966e2*log(T) ) );
		double lnPdissw = -1.6444866e3
						- 0.1374178*T
						+ 5.4979866e4/T
						+ 2.64118188e2*std::log(T)
						+ (S*1000.*1.005)*( 1.1178266e4
											+ 7.67420344*T
											- 4.515213e-3*T*T
											- 2.04872879e5/T
											- 2.17246046e3*std::log(T)
											)
		   	   	   	    + (S*1000.*1.005)*(S*1000.*1.005)*(1.70484431e2
		   	   	   	    					  	  	  	  + 0.118594073*T
														  - 7.0581304e-5*T*T
														  - 3.09796169e3/T
														  - 33.2031996*std::log(T)
		   	   	   	    					  	  	  	  ) ;
//		double lnPdissw = -1208.1964 + 44097.0/T + 186.7594*std::log(T);
		double P_eq = 1.e6*std::exp(lnPdissw); //[Pa]

		return P_eq/characteristicValue.P_c; /*ndim*/
	}

	//rate constant for hydrate dissociation;
	double DissociationRateConstant (double T) const {

		double kd_diss = parameter.HydrateDissociationRateConstant(); //defined in mol/m².Pa.s
		return kd_diss;
	}

	//rate constant for hydrate reformation;
	double FormationRateConstant_ingas (double T) const {

		double kd_form = parameter.HydrateFormationRateConstant(); //defined in mol/m².Pa.s
		return kd_form;

	}

	double FormationLimitFactor( double Sh , double Sw, double porosity ) const {

		double term = Sw*(1.-Sw-Sh) ;
		return term;

	}

	double DissociationLimitFactor( double Sh , double Sw, double porosity ) const {

		double term = Sh;
		return term  ;
	}

	//specific reaction area of hydrate in the sediment:
	double SpecificSurfaceArea (double Sh, double porosity, double permeability ) const {

		double A_s;

 		double M_base = 1.e5;
 		double M_SF = pow( porosity*(1.-Sh), 3./2. );
 		A_s = M_base * M_SF ;

// 		std::cout<< "A_s = " << A_s << std::endl ;

		return A_s;
	}

	// rate of gas generation:
	double GasGenerationRate ( double T/*K*/, double Pg/*Pa*/, double Sh,  double Sw, double XCH4,
							   double zCH4, double S, double porosity, double permeability/*m^2*/ ) const {

		double gas_gen = 0.0;

		double Peq/*Pa*/ = EquilibriumPressure( T,S )*characteristicValue.P_c;

		double potential_P = Peq/Pg - 1.  ;

		if(potential_P > 0. ){
			gas_gen =   DissociationRateConstant( T )
					  * methane.MolarMass()
					  * SpecificSurfaceArea( Sh, porosity, permeability )
					  * DissociationLimitFactor( Sh, Sw, porosity )
					  * Pg
					  * potential_P
					  ;
		}
		else if(potential_P < 0. ){
			gas_gen =   FormationRateConstant_ingas( T )
					  * methane.MolarMass()
					  * SpecificSurfaceArea( Sh, porosity, permeability )
					  * FormationLimitFactor( Sh, Sw, porosity )
					  * Pg
					  * potential_P
					  ;
		}

	    return gas_gen; /*[kg/m³s]*/
	}

	// rate of water generation:
	double WaterGenerationRate ( double gasGenRate /*kg.m^-3.s^-1*/ ) const {
      double water_gen =  gasGenRate * hydrate.HydrationNumber() * ( water.MolarMass() / methane.MolarMass() ) ;
      return water_gen;	/*[kg/m³s]*/
	}

	// rate of hydrate dissociation:
	double HydrateDissociationRate( double gasGenRate/*kg.m^-3.s^-1*/ ) const {
      double hyd_decomp= - gasGenRate * ( hydrate.MolarMass() / methane.MolarMass() ) ;
      return hyd_decomp;/*[kg/m³s]*/
	}

	// heat of hydrate dissociation:
	double HeatOfDissociation( double gasGenRate/*kg.m^-3.s^-1*/, double T/*K*/ ) const {
      double Q_decomp/*[W/m³]*/= - ( gasGenRate / hydrate.MolarMass() )
      						     * ( 56599.0 - 16.744*( T ) )
								 * 1.;

      return Q_decomp;/*[W/m³]*/
	}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
