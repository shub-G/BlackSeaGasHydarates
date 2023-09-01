class Methane
{
private:
	CharacteristicValues characteristicValue;

public:

	/* http://en.wikipedia.org/wiki/Gas_constant */
	constexpr static double Ru = 8.314462175; /* [J*mol^-1*K^-1] */

	double MolarMass() const {
		return 16.04 * 1.0e-3; 	/* [kg/mol] */
	}

	double AccentricityFactor()const {
		return 0.011 	 ;
	}

	double CriticalTemperature( ) const {
		return -82.7 + 273.15 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 45.96 * 1.0e5 ; /* [Pa] */
	}

	double Density(double T/*K*/, double Pg/*Pa*/, double z_CH4) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		double R_CH4 = Ru/MolarMass();
		rho = Pg/( z_CH4 * R_CH4 * T);

		return rho/characteristicValue.density_c; /*ndim*/
	}

	double DynamicViscosity(double T/*K*/, double Pg/*Pa*/) const {

		double mu;
		/* mu: unit -> Pa.s */

		// Sutherland Correlation:
		// ref: http://portal.tpu.ru/SHARED/n/NATASHA/Material/Tab3/Glava_1.pdf
		double C = 162; // empirical constant
		double mu_0 = 1.0707e-5;
		// ref for mu_0 :http://www.pipeflowcalculations.com/tables/gas.php
		mu_0 *= (  1. - (1./(1.0707e-5)) * (  4.8134e-14 * Pg
											- 4.1719e-20 * Pg * Pg
											+ 7.3232e-28 * Pg * Pg * Pg )
				); // Pa.s -> ref: Frauenhofer Comsol Model
		mu = mu_0 * (273.15 + C ) * ( pow( (T/273.15), 1.5) / ( T + C ) ) ;

		return mu/characteristicValue.viscosity_c;/*ndim*/

	}

	double ThermalConductivity( double T/*K*/, double Pg/*Pa*/) const {

		double kth; /* [W*m^-1*K^-1] */

#ifdef FRAUENHOFER_MODEL
		// REFERENCE: Frauenhofer Comsol Model
		kth = 43.82e-3 ; /* [W*m^-1*K^-1] */
#else
		// REFERENCE: " Thermal Conductivity of Methane for temperatures between 110 K and 310 K with Pressures upto 70 MPa
		//			  author : H.M. Roder
		//			  Journal: Journal of Thermophysics, volume 6, No 2, pages 119-142
		// assumption: dilute gas , therefore, density effects are neglected
		double A0 = - 0.8863333440 * 1.e-2 ;
		double A1 =   0.2419639784 * 1.e-3 ;
		double A2 = - 0.6997019196 * 1.e-6 ;
		double A3 =   0.1224609018 * 1.e-8 ;
		kth = A0 + A1 * T + A2 * T*T + A3 * T*T*T ;
#endif

		return kth/characteristicValue.thermalconductivity_c;/*ndim*/
	}

	double Cp_ideal( double T/*K*/, double Pg/*Pa*/ ) const {

		double Cp_i;
		/* [J/(kg*K)] */

#ifdef FRAUENHOFER_MODEL
		// REFERENCE: Frauenhofer COMSOL model
		Cp_i = 3274.0; /* [J/(kg*K)] */
#else
		/* REF: 1D Modelling of Hydrate Decomposition in Porous Media, by F. Esmailzadeh, M.E. Zeighami, J. Fathi */
		double A = 1.238 ;
		double B = 0.00313 ;
		double C = 7.905*1.0e-7 ;
		double D = -6.858*1.0e-10 ;
		Cp_i = ( A + B*T + C*T*T +D*T*T*T ) * 1000.0 ; /* [J/(kg*K)] */
#endif

		return Cp_i ; /* [J/(kg*K)] */
	}

	double Cp_res( double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {

		double Cp_res;
		/* [J/(kg*K)] */

#ifdef FRAUENHOFER_MODEL
		Cp_res = 0.;
#else
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = AccentricityFactor();
		double Tc	 = CriticalTemperature();
		double Pc 	 = CriticalPressure();

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double b = 0.07780 * ( Ru * Tc / Pc );
		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;

		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));

		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
		double B = ( b * Pg ) / ( Ru * T ) ;
		double M = ( z_CH4*z_CH4 + 2*B*z_CH4 - B*B ) / ( z_CH4 - B ) ;
		double N = ( da_T * B ) / ( b * Ru );

		Cp_res =   dda_T * ( T / (2*sqrt(2) * b ) ) * log((z_CH4+(sqrt(2)+1)*B)/(z_CH4-(sqrt(2)-1)*B))
				+ ( Ru * pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) )
				- Ru ;
#endif

		return Cp_res; /*[J/(kg.K)]*/
	}

	double Cp( double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double Cp;
		/* [J/(kg*K)] */

		Cp = Cp_ideal( T, Pg ) + Cp_res( T, Pg, z_CH4 );		/* [J/(kg*K)] */

		return Cp/characteristicValue.specificheat_c; /*ndim*/
	}

	double Cv(double T/*K*/, double Pg/*Pa*/, double z_CH4 ) const {

		double Cv;
		/* [J/(kg*K)] */

		Cv/*[J/(kg.K)]*/ =   Cp( T, Pg, z_CH4 )*characteristicValue.specificheat_c; // NOTE: Cases are checked in Cp_ideal function!

#ifdef FRAUENHOFER_MODEL
		Cv += (-1.) * Ru/MolarMass();		/* [J/(kg*K)] */
#else
		// Based on Peng Robinson's EoS
		// REFERENCE:
		double omega = AccentricityFactor();
		double Tc	 = CriticalTemperature();
		double Pc 	 = CriticalPressure();

		double kappa = 0.;
		if( omega <= 0.49 ){
			kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega ;
		}
		else{
			kappa = 0.379642 + 1.48503 * omega - 0.164423 * omega * omega ;
		}

		double ac = pow( (1 + kappa * ( 1 - sqrt(T/Tc) ) ) , 2 );

		double b = 0.07780 * ( Ru * Tc / Pc );
		double a_T = 0.45724 * ( ( Ru * Ru * Tc * Tc ) / Pc ) * ac ;

		double da_T = kappa * ac * ( ( kappa / Tc ) - (( 1 + kappa )/( sqrt( T*Tc ) ) ) );
		double dda_T = ( kappa * ac * ( 1. + kappa ) ) / ( 2. * sqrt( T*Tc ));

		double A = ( a_T * Pg ) / ( pow( ( Ru * T ) , 2 ) ) ;
		double B = ( b * Pg ) / ( Ru * T ) ;
		double M = ( z_CH4*z_CH4 + 2*B*z_CH4 - B*B ) / ( z_CH4 - B ) ;
		double N = ( da_T * B ) / ( b * Ru );

		double nonidealfactor = ( pow( M-N , 2 ) ) / ( M*M - 2.*A * (z_CH4+B) );

		Cv += (-1.) * Ru * nonidealfactor ;
#endif

		return Cv/characteristicValue.specificheat_c;/*ndim*/
	}

	double EquilibriumXCH4( double T/*K*/, double Pg/*Pa*/, double z_CH4, double S ) const {

		// REF: SUGAR TOOLBOX
		Water water;
		Salt salt;
		double Tc /*[K]*/   = water.CriticalTemperature();	//critical temperature of water
		double Pc /*[MPa]*/ = water.CriticalPressure()/1e6;	//critical pressure of water
		double tau = 1 - T/Tc;    							//dimensionless temperature
		double Ts = T/Tc;         							//reduced temperature

		// vapor pressure of water
		double a1 = -7.85951783;
		double a2 =  1.84408259;
		double a3 = -11.7866497;
		double a4 =  22.6807411;
		double a5 = -15.9618719;
		double a6 =  1.80122502;

		double xc = S * (water.MolarMass()/salt.MolarMass()) ;
		double u = 0.024;
		double c1 = 0.128746;
		double c2 = -0.731097;
		double c3 = -315.058;
		double b1 = 392.767;
		double b2 = -2464.4;
		double h=0.;
		if (xc < u){
			h = c2*xc/(xc+c1*c1) + c3*xc*xc;
		}else{
			h = ((c1*c1*c2)/std::pow(u+c1*c1 , 2.0) + 2.*c3*u)*(xc-u) + b1*(xc- u)*(xc-u)
					+ b2*(xc-u)*(xc*xc-u*u) + c2*u/(u+c1*c1) + c3*u*u;
		}

		double lnppc =  Tc/T * (  a1*tau
								+ a2*std::pow(tau,1.5)
								+ a3*std::pow(tau,3)
								+ a4*std::pow(tau,3.5)
								+ a5*std::pow(tau,4)
								+ a6*std::pow(tau,7.5)
								);
		lnppc = lnppc + h;
		double pvap /*[MPa]*/ = Pc * exp(lnppc);
		pvap /*bar*/ *= 10.;

		//parametric calculation of the chemical potential
		double d1 = 8.3143711;
		double d2 = -0.00072772168;
		double d3 = 2148.9858;
		double d4 = -0.000014019672;
		double d5 = -667434.49;
		double d6 = 0.0076985890;
		double d7 = -0.0000050253331;
		double d8 = -3.0092013;
		double d9 = 484.68502;
		double P/*bar*/ = Pg*1.e-5;
		double uch4 = d1 + d2*T + d3/T + d4*T*T + d5/(T*T) + d6*P + d7*P*T + d8*P/T + d9*P/(T*T);

		//calculation of water molar volume, approximated as saturated liquid phase volume in cm^3/mol
		double mmolh2o = 18.0153; //in g/mol
		double rhoc = 322.0; //in kg/m^3 or g/dm^3 or mg/cm^3 or ...
		double f1 = 1.99274064;
		double f2 = 1.09965342;
		double f3 = -0.510839303;
		double f4 = -1.75493479;
		double f5 = -45.5170352;
		double f6 = -674694.45;
		double rhoh2o = rhoc*(1 + f1*std::pow(tau,(1./3.))
								+ f2*std::pow(tau,(2./3.))
								+ f3*std::pow(tau,(5./3.))
								+ f4*std::pow(tau,(16./3.))
								+ f5*std::pow(tau,(43./3.))
								+ f6*std::pow(tau,(110./3.))
							  );
		double mvh2o = (mmolh2o/rhoh2o)*1000.;

		//calculation of the fugacity coefficient of water
		double g1 = -0.0142006707;
		double g2 =  0.0108369910;
		double g3 = -0.00000159213160;
		double g4 = -0.0000110804676;
		double g5 = -3.14287155;
		double g6 = 0.00106338095;
		double fch2o = std::exp(g1 + g2*P + g3*P*P + g4*P*T + g5*P/T + g6*P*P/T);

		//
		double Tr = T/190.6;
		double Pr = Pg/(46.41e5);
		double Vr = z_CH4*Tr/Pr;
		double k1=0.0872553928,
			   k2=-0.752599476,
			   k3=0.375419887,
			   k4=0.0107291342,
			   k5=0.00549626360,
			   k6=-0.0184772802,
			   k7=3.18993183e-4,
			   k8=2.11079375e-4,
			   k9=2.01682801e-5,
			   k10=-1.65606189e-5,
			   k11=1.19614546e-4,
			   k12=-1.08087289e-4,
			   Alpha=0.0448262295,
			   Beta=0.75397,
			   Gamma=0.077167;
		double B = k1 + k2/(Tr*Tr) + k3/(Tr*Tr*Tr);
		double C = k4 + k5/(Tr*Tr) + k6/(Tr*Tr*Tr);
		double D = k7 + k8/(Tr*Tr) + k9/(Tr*Tr*Tr);
		double E = k10 + k11/(Tr*Tr) + k12/(Tr*Tr*Tr);
		double F = Alpha/(Tr*Tr*Tr);
		double G = F/(2.*Gamma)*( Beta + 1.0 - (Beta + 1.0 + (Gamma/Vr*Vr) ) * std::exp(-Gamma/(Vr*Vr)) );
		double lnPHI = z_CH4 - 1.0 - std::log(z_CH4) + B/Vr + C/(2.0*Vr*Vr) + D/(4.0*Vr*Vr*Vr*Vr) + E/(5.0*Vr*Vr*Vr*Vr*Vr) + G;
		double PHI = std::exp(lnPHI);

		// Calculation of molality of CH4 in solution
		double R2006 = 83.14472; //gas constant in bar cm^3/mol/K for Duan 2006
		double mfh2o = 1.0; //mfh2o = mh2o./(mh2o + mNa + mCl + mK + mCa + mMg + mSO4)
		//---> calculation of mole fraction of water in gas phase
		double ysw = mfh2o * (pvap / fch2o / P) * std::exp( mvh2o*(P-pvap)/R2006/T );
		//---> calculation of mole fraction of CH4 in gas phase
		double ych4sw = 1. - ysw;
		//--->
		double lnmch4 = std::log(ych4sw*PHI*P) - uch4;
		double mCH4sw = std::exp(lnmch4); //Molality of dissolved CH4 [mol/kg]
		double xch4sw = mCH4sw * water.MolarMass(); //mol/kg * kg/mol
		//--->

		//Calculation of dissociation pressure of methane hydrate in seawater
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
		double pdissw = std::exp(lnPdissw); //[MPa]

		//Calculation of solubilities in seawater
		//---> at the hydrate dissociation pressure Pdissw (-> for 3 coexisting phases: gas, hydrate, seawater)
		double lnch4ghswPdis = -2.5640213e5
							 - 1.6448053e2*T
							 + 9.1089042e-2*T*T
							 + 4.90352929e6/T
							 + 4.93009113e4*std::log(T)
							 + S*1000.* (-5.16285134e2
									     - 0.33622376*T
										 + 1.88199047e-4*T*T
										 + 9.76525718e3/T
										 + 9.9523354e1*std::log(T)
							 	 	 	 );
		//---> at hydrostatic pressure (-> seawater in equilibrium with methane hydrate; no gas)
		double lnch4ghsw = lnch4ghswPdis
						 + (Pg*1.e-6-pdissw)*(5.04597e-2
								 	 	 	 + 7.64415e-4*S*1000.
											 - T*(3.90236e-4 + 5.48947e-6*S*1000.0)
											 + T*T*(7.06154e-7 + 9.87742e-9*S*1000.0)
											 )
						+ (Pg*1.e-6-pdissw)*(Pg*1.e-6-pdissw)*(   7.57285e-5
																- 1.90867e-8*S*1000.
																- 1.4483e-10*S*S*1000.*1000.
																- T*(1.96207e-7 - 6.67456e-11*S*1000.0)
															   );
		double mCH4ghsw = std::exp(lnch4ghsw); //[mol/kg]
		double xch4ghsw = mCH4ghsw * water.MolarMass(); //mol/kg * kg/mol

		double xcheqb = 0.;
	    if(Pg*1.e-6 > pdissw) xcheqb=xch4ghsw;
	    else xcheqb = 0.75*xch4sw;

		return xcheqb; /*ndim*/

	}

};
