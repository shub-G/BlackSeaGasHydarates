template<typename GV, typename Parameters, typename PGMap>
class Soil
{
private:
	  const GV& gv;
	  const Parameters& parameter;
	  const PGMap& pgmap ;
	  
	  CharacteristicValues characteristicValue;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  Soil ( const GV& gv_ , const Parameters& parameter_ , const PGMap& pgmap_ )
  : gv( gv_ ), parameter(parameter_), pgmap( pgmap_ )
  {}


	double SedimentPorosity
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.layer_ztop();
		auto prop_L = parameter.layer_properties();

		double por = 0.;

		if( x[dim-1]>z_L[0] ){
			por = prop_L[0][0];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i] and x[dim-1]<z_L[i-1] ){
				por = prop_L[i][0];
			}
		}

		return por;
	}

	double SedimentPermeability
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);
		
		auto z_L = parameter.layer_ztop();
		auto prop_L = parameter.layer_properties();

		double K = 0.; /*m^2*/

		if( x[dim-1]>z_L[0] ){
			K = prop_L[0][1];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i] and x[dim-1]<z_L[i-1] ){
				K = prop_L[i][1];
			}
		}

		return K/characteristicValue.permeability_c; /*ndim*/
	}

	// vector coefficient
	Dune::FieldVector<double,dim>
	SedimentPermeabilityVector
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		Dune::FieldVector<double,dim> x = element.geometry().global(xlocal);

		auto z_L = parameter.layer_ztop();
		auto prop_L = parameter.layer_properties();

		double Kxx = 0.; /*m^2*/
		double Kzz = 0.; /*m^2*/

		if( x[dim-1]>z_L[0] ){
			Kxx = prop_L[0][1];
			Kzz = prop_L[0][8];
		}
		for( int i=1; i<z_L.size(); i++ ){
			if( x[dim-1]>z_L[i] and x[dim-1]<z_L[i-1] ){
				Kxx = prop_L[i][1];
				Kzz = prop_L[i][8];
			}
		}
		Dune::FieldVector<double,dim> PermeabilityVector;
		PermeabilityVector[0] = Kxx/characteristicValue.permeability_c; /*ndim*/
		PermeabilityVector[1] = Kzz/characteristicValue.permeability_c; /*ndim*/

		return PermeabilityVector; /*ndim*/
	}

	double SoilGrainRadius
	(const typename GV::Traits::template Codim<0>::Entity& element,
	 const Dune::FieldVector<double,dim>& xlocal) const {

		auto x = element.geometry().global(xlocal);

		// Bear et at, 1972
		double por  = SedimentPorosity( element,xlocal );
		double perm = SedimentPermeability( element,xlocal )*characteristicValue.permeability_c;
		double rp = sqrt( 45.0 * perm * pow( 1- por , 2.0 )/pow( por,3.0) );
		return rp/characteristicValue.x_c; /*ndim*/
	}

	double Density() const {
		/* unit -> kg/m^3 */
		double rho = 2600.0;
		return rho/characteristicValue.density_c; /*ndim*/
	}

	double ThermalConductivity() const {
		/* unit -> W/mK */
		double kth = 3.0;
		return kth/characteristicValue.thermalconductivity_c; /*ndim*/
	}

	double Cp() const {
		/* unit -> J/kg.K */
		double Cp = 1000.0;
		return Cp/characteristicValue.specificheat_c; /*ndim*/
	}

	double Cv() const {
		/* unit -> W/kg.K */
		double Cv = Cp(); /*ndim*/
		return Cv;/*ndim*/
	}

	double Tortuosity( double porosity ) const {
		return porosity * porosity ;/*ndim*/
	}

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}

};
