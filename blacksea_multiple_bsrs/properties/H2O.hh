class Water
{
private:
	CharacteristicValues characteristicValue;
	Salt salt;

public:

	double CriticalTemperature( ) const {
		return 647.096 ; /* [K] */
	}

	double CriticalPressure( ) const {
		return 22.064 * 1.0e6 ; /* [Pa] */
	}

	double MolarMass( ) const {
		/* unit -> kg/mol */
		return 18.0/1000;
	}

	double Density( double T/*K*/, double Pw/*Pa*/, double S ) const {

		double rho;
		/* rho: unit -> kg/m^3 */

		/*
		 * averages values & expansion coefficients: ρ0=1027 kg/m^3,  T0=10°C,  S_0=35 g/kg
		 * Thermal expansion: \alpha_T=0.15 kg/(m^3 °C)
		 * Salinity contraction: \alpha_S=0.78 kg/(m^3 g/kg)
		 * Pressure compressibility: \alpha_P=0.0045 kg/(m^3 dbar)
		 * UNESCO EOS-80 : Equation of state for seawater
		 * We use a linear EOS (web ref:http://www.ccpo.odu.edu/~atkinson/OEAS405/Chapter2_Stratified_Ocean/Lec_04_DensityEOS.pdf)
		 */

		double rho_0 = 1027.0;
		double T_0 = 10.;
		double S_0 = 0.035;
		double alpha_T = -0.15;
		double alpha_S = 0.78*1e3;
		double alpha_P = 0.0045;

		rho = rho_0
			+ (   alpha_P*(Pw*1.e-4)
				+ alpha_T*((T-273.15)-T_0)
				+ alpha_S*(S-S_0)
			  );

		return rho/characteristicValue.density_c; /*ndim*/
	}

	double DynamicViscosity( double T/*K*/, double Pw/*Pa*/, double S ) const {
		double mu;
		/* mu: unit -> Pa.s */

#ifdef FRAUENHOFER_MODEL
		mu = 0.5 * 1.0e-3 ;
#else
		// REFERENCE:
		double mu_0 = 0.001792 ; // kg/m/s
		double a = - 1.94 ;
		double b = - 4.80 ;
		double c =  6.74 ;
		double T0 = 273.15 ; // K
		double Tr = T0/T;
		mu = mu_0 * exp( a + b * Tr + c * Tr*Tr );
#endif

		return mu/characteristicValue.viscosity_c;/*ndim*/
	}

	double ThermalConductivity( double T/*K*/, double Pw/*Pa*/, double S ) const {

		double kth;
		/* kth: unit -> W.m^-1 K^-1 */

		kth = 0.57153*( 1 + 0.003*(T-273.15) - 1.025e-5*(T-273.15)*(T-273.15) + 6.53e-10*Pw - 0.29*S );

		return kth/characteristicValue.thermalconductivity_c; /*ndim*/
	}

	double Cp( double T/*K*/, double Pw/*Pa*/, double S ) const {
		double Cp;
		/* Cp: unit -> J*kg^-1*K^-1 */

		Cp = 3945.0 ;

		return Cp/characteristicValue.specificheat_c;/*ndim*/

	}

	double Cv( double T/*K*/, double Pw/*Pa*/, double S ) const {
		double Cv;
		/* mu: unit -> J*kg^-1*K^-1 */

		Cv = Cp( T, Pw, S ); /*ndim*/

		return Cv; /*ndim*/

	}

	double EquilibriumYH2O( double T/*K*/, double Pg/*Pa*/, double S ) const {

		double psat;   /* [Pa] */

		// REF: SUGAR TOOLBOX

		double Tc /*[K]*/  = CriticalTemperature();	//critical temperature of water
		double Pc /*[Pa]*/ = CriticalPressure();	//critical pressure of water
		double tau = 1 - T/Tc;    					//dimensionless temperature
		double Ts = T/Tc;         					//reduced temperature

		// vapor pressure of water
		double a1 = -7.85951783;
		double a2 =  1.84408259;
		double a3 = -11.7866497;
		double a4 =  22.6807411;
		double a5 = -15.9618719;
		double a6 =  1.80122502;

		double xc = S * (MolarMass()/salt.MolarMass()) ;
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

		// Calculation of molality of CH4 in solution
		double R2006 = 83.14472; //gas constant in bar cm^3/mol/K for Duan 2006
		double mfh2o = 1.0; //mfh2o = mh2o./(mh2o + mNa + mCl + mK + mCa + mMg + mSO4)
		//---> calculation of mole fraction of water in gas phase
		double ysw = mfh2o * (pvap / fch2o / P) * std::exp( mvh2o*(P-pvap)/R2006/T );

		return ysw; /*ndim*/
	}

};
