template<typename GV, typename Parameters, typename PGMap>
class HydraulicProperties
{
private:
	  const GV& gv;
	  const Parameters& parameter;
	  const PGMap& pgmap ;
	  
	  const static int dim = GV::dimension;

	  const static int numOfParams  = 6;
	  const static int id_Pentry 	= 0;
	  const static int id_lambda 	= 1;
	  const static int id_Swr 		= 2;
	  const static int id_Sgr 		= 3;
	  const static int id_m			= 4;
	  const static int id_beta		= 5;

	  Soil<GV,Parameters,PGMap> soil;
	  CharacteristicValues characteristicValue;

public:

  //! construct from grid view
  HydraulicProperties ( const GV& gv_ , const Parameters& parameter_ , const PGMap& pgmap_ )
: gv( gv_ ),
  parameter(parameter_),
  pgmap(pgmap_),
  soil(gv_,parameter_,pgmap_)
{}

	/* PARAMETERS FOR THE HYDRAULIC PROPERTY CONSTITUTIVE LAW (BROOKS-COREY) */
	std::vector<double>
	BrooksCoreyParameters( const typename GV::Traits::template Codim<0>::Entity& element,
	   	   	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		std::vector<double> BCParams (numOfParams,0.);

		auto z_L = parameter.layer_ztop();
		auto prop_L = parameter.layer_properties();

		if( x[dim-1]>z_L[0] ){
			BCParams[id_Pentry] = prop_L[0][2] ; /*Pa*/
			BCParams[id_lambda] = prop_L[0][3] ;
			BCParams[id_Swr] 	= prop_L[0][4] ;
			BCParams[id_Sgr] 	= prop_L[0][5] ;
			BCParams[id_m] 		= prop_L[0][6] ;
			BCParams[id_beta] 	= prop_L[0][7] ;
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i] and x[dim-1]<z_L[i-1] ){
				BCParams[id_Pentry] = prop_L[i][2] ; /*Pa*/
				BCParams[id_lambda] = prop_L[i][3] ;
				BCParams[id_Swr] 	= prop_L[i][4] ;
				BCParams[id_Sgr] 	= prop_L[i][5] ;
				BCParams[id_m] 		= prop_L[i][6] ;
				BCParams[id_beta] 	= prop_L[i][7] ;
			}
		}

		return BCParams;
	}

	/* EFFECTIVE SATURATION */

	double EffectiveSw( double Sw,
			 	 	 	double Sh,
						double Swr,
						double Sgr ) const {


		double Sw_max = 1. - Sh - Sgr;
		double Swe = ( Sw - Swr )/( Sw_max - Swr );

		if( Swe>1.){
			Swe=1.;
		}
		if( Swe<0.){
			Swe=0.;
		}

		return Swe;
	}

	double dSwe_dSw( double Sw,
			 	 	 double Sh,
					 double Swr,
					 double Sgr) const {

		double Sw_max = 1. - Sh - Sgr;
		double dSwe =  1./( Sw_max - Swr );

		return dSwe;
	}

	double dSwe_dSh( double Sw,
					 double Sh,
					 double Swr,
					 double Sgr) const {

		double Sw_max = 1. - Sh - Sgr;
		double Swe = ( Sw - Swr )/( Sw_max - Swr );
		double dSwe = Swe /( Sw_max - Swr );

		return dSwe;
	}

	/* SUCTION/CAPILLARY PRESSURE */

	double CapillaryPressure(  const typename GV::Traits::template Codim<0>::Entity& element,
  	   	   	 	 	 	 	   const Dune::FieldVector<double,dim>& xlocal ,
							   double Sw,
							   double Sh,
							   double porosity) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);

		double Pentry/*Pa*/	= BCParams[id_Pentry];
		double lambda 		= BCParams[id_lambda];
		double Sgr 			= BCParams[id_Sgr];
		double Swr 			= BCParams[id_Swr];
		double m 			= BCParams[id_m];
		double beta 		= BCParams[id_beta];

		double eta = (1/lambda);
		double Swe = EffectiveSw( Sw,Sh,Swr,Sgr );
		double Pc = 0.; /*Pa*/
		double a = 0.05 ;

		if( Swe > a ){
			Pc = Pentry * pow( Swe, -eta );
		}
		else if ( Swe <= a ){
			double Pc_a /*Pa*/  = Pentry * pow( a, -eta );
			double dPc_a /*Pa*/ = dPc_dSwe( a,Pentry,lambda ) ;
			Pc/*Pa*/ = Pc_a/*Pa*/ + dPc_a/*Pa*/ * ( Swe - a );
		}
		else {
			std::cout<< " ERROR in " << __FILE__
					 << " function: SuctionPressure( element,xlocal,Sw,Sh )"
					 << "  , Swe = " << Swe
					 << "  , Sw  = " << Sw
					 << "  , Sh  = " << Sh
					 << "  , Pc  = " << Pc << std::endl;
//			exit(0);
		}

		if( Pc < -1e-3 ){
			std::cout<< " Pc is -ve " << std::endl;
			std::cout<< " ERROR in " << __FILE__
					 << " function: SuctionPressure( element,xlocal,Sw,Sh )"
					 << "  , Swe = " << Swe
					 << "  , Sw  = " << Sw
					 << "  , Sh  = " << Sh
					 << "  , Pc  = " << Pc << std::endl;
//			exit(0);
		}

		double SF_Sh = PcSF1( Sh,lambda,m );

		double porosity_0 = soil.SedimentPorosity( element,xlocal );
		double SF_por = PcSF2( porosity, porosity_0, beta );

		Pc *=SF_Sh*SF_por; /*Pa*/
		return Pc/characteristicValue.P_c; /*ndim*/
	}

	double dPc_dSwe( double Swe,
					 double Pentry, /*Pa*/
					 double lambda ) const {

		double eta = (1/lambda);
		double dPc = 0.; /*Pa*/
		double a = 0.05 ;

		if( Swe > a ){
			dPc/*Pa*/ = Pentry * (-1./lambda) * std::pow( Swe , -(1./lambda) - 1. );
		}
		else if ( Swe <= a ){
			double dPc_a  = Pentry * (-1./lambda) * std::pow( a , -(1./lambda) - 1. ) ;
			double ddPc_a = Pentry * (-1./lambda) * (-1./lambda-1.) * std::pow( a , -(1./lambda) - 2. );
			dPc/*Pa*/ = dPc_a + ddPc_a * ( Swe - a );
		}
		else {
			std::cout<< " ERROR in HydraulicProperties::dPc_dSwe( Swe ) "
					 << "  , Swe = " << Swe
					 << "  , dPc  = " << dPc << std::endl;
//			exit(0);
		}

		return dPc; /*Pa*/
	}

	/* Pc SCALING FACTORS */

	double PcSF1( double Sh,
				  double lambda,
				  double m) const {

		double PcSF=0. ;
		double eta = std::abs((lambda*m - 1.)/(lambda*m));
		double a = 0.95;

		if( Sh < a ){
			PcSF =  pow( 1.0 - Sh , -eta );
		}
		else if( Sh >= a ){
			double PcSF_a = std::pow( 1.0 - a , -eta );
			double dPcSF_a = eta * std::pow( 1.0 - a , -eta-1 );
			PcSF = PcSF_a + dPcSF_a * ( Sh - a );
		}
		else {
			std::cout<< " ERROR in " << __FILE__
					 << " function: PcSF1( e, xlocal, Sh ). "
					 << " Sh = " << Sh
					 << std::endl;
//			exit(0);
		}
		return PcSF ;
	}

	double PcSF2( double phi,
				  double phi_0,
				  double beta) const {

		double PcSF;

		double a = 0.05 ;
		if ( phi > a ){
			double term = ( phi_0/phi ) * ( ( 1-phi )/( 1. - phi_0 ));
			PcSF =pow( term, beta );
		}
		else if( phi <= a ){
			double term_a = ( phi_0/a ) * ( ( 1-a )/( 1. - phi_0 ));
			double PcSF_a = std::pow( term_a,beta );
			double dPcSF_a = 0.;
			double C = std::pow( phi_0/( 1. - phi_0 ) , beta );
			dPcSF_a -= beta * C * std::pow( 1.-a , beta-1. ) * std::pow( a , -beta-1. );
			PcSF = PcSF_a + dPcSF_a * ( phi - a );
		}
		else {
			std::cout<< " ERROR in " << __FILE__
					 << " function: PcSF2( e, xlocal, phi_0, phi )" << std::endl;
//			exit(0);
		}

		return PcSF ;
	}

	/* RELATIVE PERMEABILITIES */

	double krw( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
				double Sw,
				double Sh ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];

		double Swe = EffectiveSw(Sw,Sh,Swr,Sgr);

		double kr = std::pow(Swe, (2.0/lambda + 3.0) );
		if( Swe>1.){
			kr=1.;
		}
		if( Swe<0.){
			kr=0.;
		}

		return kr ;
	}

	double krg( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	const Dune::FieldVector<double,dim>& xlocal ,
				double Sw,
				double Sh ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];

		double Swe = EffectiveSw(Sw,Sh,Swr,Sgr);

		double kr = std::pow(1.0-Swe, 2.0) * ( 1.0 - std::pow(Swe, (2.0/lambda + 1.0) ) );

		if( Swe>1.){
			kr=0.;
		}
		if( Swe<0.){
			kr=1.;
		}
		return kr;
	}

	/* PERMEABILITY SCALING FACTORS */

	double PermeabilityScalingFactor( const typename GV::Traits::template Codim<0>::Entity& element,
	   	 	 	 	 	 	 	 	  const Dune::FieldVector<double,dim>& xlocal ,
									  double Sh,
									  double porosity ) const {

		auto BCParams = BrooksCoreyParameters(element,xlocal);
		double Pentry 	= BCParams[id_Pentry];
		double lambda 	= BCParams[id_lambda];
		double Sgr 		= BCParams[id_Sgr];
		double Swr 		= BCParams[id_Swr];
		double m 		= BCParams[id_m];
		double beta 	= BCParams[id_beta];

		double SF_Sh = KSF1( Sh,m );

		double porosity_0 = soil.SedimentPorosity( element,xlocal );
		double SF_por = KSF2( porosity, porosity_0, beta );

		double SF =SF_Sh*SF_por;
		return SF;
	}

	double KSF1( double Sh,
				 double m ) const {

		return std::pow( (1.0 - Sh) , (5*m + 4)/(2*m) );
	}

	double KSF2( double phi,
				 double phi_0,
				 double beta ) const {
		// POWER LAW MODEL PROPOSED BY CIVAN (2001)
		// Read " kinetic simulation of methane hydrate formation and issociation in porous media " by Xuefei Sun and Kishore Mohanty
		// phi_0 = basePorosity_initial (and NOT sediemntPorosity!! )
		// phi is basePorosity at current time

		double term1 = phi / phi_0 ;
		double term2 = ( 1-phi_0 ) / ( 1. - phi ) ;

		double KSF=0.;
		double a = 0.95 ;

		if( phi < a ){
			KSF = term1 * std::pow( ( term1 * term2 ) , 2*beta );
		}
		else if ( phi >= a ){
			double C = std::pow( 1-phi_0 , 2*beta ) / std::pow( phi_0 , 2*beta +1 );
			double KSF_a = C * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta );
			double dKSF_a =    C * ( 2*beta +1 ) * std::pow( a, 2*beta   ) * std::pow( 1.-a , -2*beta      )
							 + C * ( -2*beta   ) * std::pow( a, 2*beta+1 ) * std::pow( 1.-a , -2*beta - 1. );
			KSF = KSF_a + dKSF_a * ( phi - a );
		}
		else {
			std::cout<< " ERROR in " << __FILE__
					 << " function: KSF2( e, xlocal, phi_0, phi )" << std::endl;
//			exit(0);
		}
		return KSF;
	}


  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
