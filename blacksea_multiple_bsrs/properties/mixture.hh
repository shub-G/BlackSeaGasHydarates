class Mixture
{
private:
	Methane methane;
	Water water;
	Salt salt;
#ifdef FRAUENHOFER_MODEL
  	FrauenhoferFunction eos;
#elif defined(PENG_ROBINSON_EOS)
  	PengRobinson eos;
#else
  	BaseEoS eos;
#endif
	CharacteristicValues X_c;

public:

	/* MOLE FRACTIONS ( X -> liquid ; Y -> gas ) */

	std::vector<double> EquilibriumMoleFractions( double T/*K*/, double Pg/*Pa*/, double Xc, double z )const{

		double Y_H2O, Y_CH4, X_H2O, X_CH4;

		// NOTE: it is not necessary to check state dependence for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/water.MolarMass());
		Y_H2O = water.EquilibriumYH2O(T,Pg,S);
		Y_CH4 = 1.-Y_H2O;
		X_CH4 = methane.EquilibriumXCH4(T,Pg,z,S);
		X_H2O = 1. - Xc - X_CH4;

		std::vector<double> X(Indices::numOfComps,0.);
		X[Indices::compId_XCH4] = X_CH4;
		X[Indices::compId_XH2O] = X_H2O;
		X[Indices::compId_YCH4] = Y_CH4;
		X[Indices::compId_YH2O] = Y_H2O;

		return X;
	}

	double YCH4( double X_CH4, double T, double Pg, double Xc, double z )const{

		double Y_CH4;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/water.MolarMass());
		Y_CH4 = 1.-water.EquilibriumYH2O(T,Pg,S);

		return Y_CH4;
	}

	double XH2O( double Y_H2O, double T, double Pg, double Xc )const{

		double X_H2O;

		// NOTE: it is not necessary to check case1,2 for fncs f_CH4 and f_H2O because the cases are already determined within classes CH4 and H2O.
		double S = Xc * (salt.MolarMass()/water.MolarMass());
		auto zCH4 = eos.EvaluateCompressibilityFactor( T,Pg );
		X_H2O = 1. - Xc - methane.EquilibriumXCH4(T,Pg,zCH4,S);

		return X_H2O;
	}

	/* MASS TRANSFER COEFFICIENTS */

	double DiffCoeffH2OInGas( double T/*K*/, double Pg/*Pa*/ ) const {

		double D;
		/* m^2/s */

		double a0 = 0. ;
		double a1 = 2.26e-9;
		double a2 = 0.002554;

		D = ( a0 + a1*T + a2/Pg ) ;

		return D/X_c.dispersivity_c; /*ndim*/

	}

	double DiffCoeffCH4InLiquid( double T/*K*/, double Pw/*Pa*/ ) const {

		double D; /* m^2/s */

		double A = 0.003475 ; 	/* K */
		double B = 1.57e-5;		/* cm^2/s */

		D = B * exp(-A/T) * 1.0e-6;

		return D/X_c.dispersivity_c; /*ndim*/
	}

};
