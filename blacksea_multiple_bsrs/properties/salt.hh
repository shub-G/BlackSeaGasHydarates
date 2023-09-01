class Salt
{
private:
	CharacteristicValues characteristicValue;

public:

	double MolarMass( ) const {
		/* unit -> kg/mol */
		return 58.4/1000;
	}

	double DiffCoeff( double T/*K*/, double Pw/*Pa*/ ) const {

		double D = 1.*1.0e-9;	/* m^2/s */
		return D/characteristicValue.dispersivity_c; /*ndim*/
	}

	double Source( ) const {
		return 0.; /*kg.m^-3.s^-1*/
	}

};
