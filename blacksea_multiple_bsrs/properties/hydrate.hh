class Hydrate {

private:
	CharacteristicValues characteristicValue;

public:

	double Density( ) const {
		/* unit -> kg/m^3 */

		double rho_h = 920.0 ;

		return rho_h/characteristicValue.density_c;/*ndim*/

	}

	double MolarMass() const
	{
		/* unit -> kg/mol */
		return 119.5/1000;
	}

	double HydrationNumber() const
	{
		return 5.90;
	}

	double ThermalConductivity( double T/*K*/, double P/*Pa*/ ) const
	{
		double kth;
		/* kth: unit -> W.m^-1 K^-1 */

		kth = 0.5 ;

		return kth/characteristicValue.thermalconductivity_c;/*ndim*/
	}

	double Cp( double T/*K*/, double P/*Pa*/ ) const
	{
		double Cp;
		/* Cp: unit -> J/kg.K */

		Cp = ( 1.9370547e-05*T*T*T - 1.5151760e-02*T*T + 3.9553876*T - 342.70565 )*1.0e3;

		return Cp/characteristicValue.specificheat_c; /*ndim*/

	}

	double Cv( double T/*K*/, double P/*Pa*/ ) const
	{
		double Cv;
		/* Cv: unit -> J/kg.K */

		Cv = Cp( T, P ) ;/*ndim*/

		return Cv;/*ndim*/

	}

};
