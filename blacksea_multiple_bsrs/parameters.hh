template<typename PTree>
class Parameters
{
private:
	const PTree& ptree;
	double *time;
	double *dt;

	MeshParameters<PTree> mesh;
	const double pi = 3.14159265358979323846;
	const static int dim = MeshParameters<PTree>::dimension;
	const double eps = 1.e-3;
	CharacteristicValues X_c;

	double Zmax;
	double Xmax;
	double HYD_thickness;
	double SHmax;
	double BSR_depth;
	std::vector<double> z;
	std::vector<std::vector<double> > prop;
	bool gravity_flag;
	double gravity_magnitude;
	double kd;
	double kf;
	double initial_salinity;
	double initial_nxch4;
	double PSF_Pw;
	double PSF_T;
	double thermal_gradient;
	
public:

  //! constructor
  Parameters (const PTree& ptree_, double *time_/*ndim*/, double *dt_/*ndim*/)
  :ptree(ptree_),
   time(time_),
   dt(dt_),
   mesh(ptree_)
  {
	/* layer numbering: example with numLayers=3
	Zmax___________
		|
		|
		|	layer0
		|
	______________ z0
		|
		|
		|	layer1
		|
	______________ z1
		|
		|
		|	layer2
		|
    z=0_______________ z2
	*/

	Zmax=ptree.get("grid.ug.LZ",(double)1000.); /*m*/
	Xmax=ptree.get("grid.ug.LX",(double)1000.); /*m*/

	HYD_thickness=ptree.get("paleo_conditions.hydrate_layer.thickness",(double)50.); /*m*/
	SHmax=ptree.get("paleo_conditions.hydrate_layer.max_sh",(double)0.30); /*-*/
	BSR_depth=ptree.get("paleo_conditions.BSR_depth",(double)360.); /*m*/
	initial_salinity = ptree.get("paleo_conditions.salinity",(double)0.035);
	initial_nxch4 = ptree.get("paleo_conditions.dissolved_methane_fraction",(double)0.95);

	int numLayers = ptree.get("sediment.number_of_layers",(int)1);
	/*z -> z-location of the bottom of a layer*/
	z = std::vector<double> (numLayers,0.);
	for(int n_layer=0; n_layer<numLayers; n_layer++ ){
		std::string name = "sediment.layer"+std::to_string(n_layer);
		z[n_layer] = ptree.get(name+".z",(double)0.); /*m*/
	}
	prop = std::vector<std::vector<double> > (numLayers,std::vector<double>(9, 0.));
	for(int n_layer=0; n_layer<numLayers; n_layer++ ){
		std::string name = "sediment.layer"+std::to_string(n_layer);
		prop[n_layer][0] = ptree.get(name+".porosity",	(double)0.5);
		prop[n_layer][1] = ptree.get(name+".permeability_xx",(double)1.e-15);/*m^2*/
		prop[n_layer][8] = ptree.get(name+".permeability_zz",(double)1.e-15);/*m^2*/
		prop[n_layer][2] = ptree.get(name+".entry_pressure",(double)1.e5); /*Pa*/
		prop[n_layer][3] = ptree.get(name+".lambda",(double)1.2);
		prop[n_layer][4] = ptree.get(name+".swr",	(double)0.);
		prop[n_layer][5] = ptree.get(name+".sgr",	(double)0.);
		prop[n_layer][6] = ptree.get(name+".m",		(double)3.);
		prop[n_layer][7] = ptree.get(name+".beta",	(double)1.);
	}

	gravity_flag = ptree.get("gravity.flag",(bool)true);
	gravity_magnitude = ptree.get("gravity.magnitude",(double)9.81);

	kd = ptree.get("hydrate_phase_change.dissociation_rate",(double)1.e-18);/*mol/m².Pa.s*/
	kf = ptree.get("hydrate_phase_change.formation_rate",(double)1.e-18);/*mol/m².Pa.s*/

	PSF_Pw = ptree.get("paleo_conditions.sea_floor_pressure",(double)10); /*MPa*/
	PSF_Pw *= 1.e6; /*convert to Pa*/

	PSF_T  = ptree.get("paleo_conditions.sea_floor_temperature",(double)4.0); /*degC*/
	PSF_T += 273.15; /*convert to K*/

	thermal_gradient = ptree.get("paleo_conditions.regional_temperature_gradient",(double)35.0); /*degC/km*/
	thermal_gradient *= 1./1000.; /*degC/m */

  }
  
	// INITIAL GAS HYDRATE LAYER
	bool isInitialHydrateLayer( Dune::FieldVector< double, dim > globalpos /*ndim*/ )const{
		double zPSFC /*ndim*/ = mesh.z_PSFC(globalpos[0]/*ndim*/);
		zPSFC *= X_c.x_c; /*m*/
		double Z_A /*m*/ = zPSFC-BSR_depth+HYD_thickness;
		double Z_B /*m*/ = zPSFC-BSR_depth;
		if( (globalpos[1]<Z_A/X_c.x_c+eps) and (globalpos[1]>Z_B/X_c.x_c-eps)){
			return true;
		}else return false;
	}

	//SEDIMENT LAYERS AND HYDRAULIC PROPERTIES
	std::vector<double> layer_ztop() const {
		return z;
	}
	std::vector< std::vector<double> > layer_properties() const {
//		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		return prop;
		/*values returned WITH dimensions*/
	}
	
	//HYDRATE REACTION KINETIC CONSTANTS
	double HydrateDissociationRateConstant() const {
//		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		return kd;
		/*mol/m².Pa.s*/
	}
	double HydrateFormationRateConstant() const {
//		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		return kf;
		/*mol/m².Pa.s*/
	}
	
	//INITIAL PALEO CONDITIONS
	double PSF_Pressure() const {
		return PSF_Pw; /*Pa*/
	}
	double PSF_Temperature() const {
		return PSF_T; /*K*/
	}
	double RegionalThermalGradient() const {
		return thermal_gradient; /*degC/m*/
	}
	double BSR_Depth() const {
		return BSR_depth; /*m*/
	}
	double PaleoSalinity() const {
		return initial_salinity; /*kg/kg*/
	}
	double InitialXCH4Fraction() const {
		return initial_nxch4; /*-*/
	}
	
	//INITIAL HYDRATE DISTRIBUTION
	double InitialHydrateSaturation(Dune::FieldVector<double,dim> globalpos /*ndim,ndim*/) const {
		  double Sh = 0.0;
		  // hydrate lens
			double zPSFC /*ndim*/ = mesh.z_PSFC(globalpos[0]/*ndim*/);
			zPSFC *= X_c.x_c; /*m*/
			double Z_A /*m*/ = zPSFC-BSR_depth+HYD_thickness;
			double Z_B /*m*/ = zPSFC-BSR_depth;
		  if( isInitialHydrateLayer(globalpos) ) {
			  double A = -Z_A*Z_B;
			  double B = Z_A+Z_B;
			  double C = -1.;
			  Sh = SHmax/std::pow(HYD_thickness/2.,2);
			  Sh *= A+B*globalpos[1]*X_c.x_c+C*globalpos[1]*X_c.x_c*globalpos[1]*X_c.x_c;
		  }

		return Sh; /*-*/
	}
	
    /* SEDIMENT BURIAL */
	double SedimentationDepth( double x /*m*/ ) const {

		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/

		double sedimentation_depth = 0.;

		double mSF1=21.43 , mSF2=80.00 , mSF3=-16.67 , mSF4=40.00 , mSF5=-9.23 ;
		double mA1=21.43  , mA2=53.33  , mA3=-36.67  , mA4=30.00  , mA5=-27.69 ;
		double mB1=21.43  , mB2=53.33  , mB3=-40.00  , mB4=40.00  , mB5=-30.77 ;
		double mC1=21.43  , mC2=53.33  , mC3=-61.67  , mC4=67.50  , mC5=-32.31 ;

		double cSF1=175.00 , cSF2=-235.00 , cSF3=586.67 , cSF4=-65.00  , cSF5=599.62 ;
		double cA1=120.00  , cA2=-103.33  , cA3=661.67  , cA4=-105.00  , cA5=673.85  ;
		double cB1= 60.00  , cB2=-163.33  , cB3=630.00  , cB4=-290.00  , cB5=665.38  ;
		double cC1= 30.00  , cC2=-193.33  , cC3=784.17  , cC4=-701.72  , cC5=646.15  ;

		double X1 = 7.0*1000. , X2=9.0*1000. , X3=12.0*1000. , X4=14.0*1000. , X5=20.0*1000. ;
		double 	t1 = 60.  * 1000. * 365.*86400. ,
				t2 = 107. * 1000. * 365.*86400. ,
				t3 = 163. * 1000. * 365.*86400. ,
				t4 = 225. * 1000. * 365.*86400. ,
				t5 = 300. * 1000. * 365.*86400. ;

		double delta_m_CB  = 0.,  delta_c_CB  = 0., delta_z_CB  = 0., delta_t_CB  = 0., vs_CB  = 0.;
		double delta_m_BA  = 0.,  delta_c_BA  = 0., delta_z_BA  = 0., delta_t_BA  = 0., vs_BA  = 0.;
		double delta_m_ASF = 0.,  delta_c_ASF = 0., delta_z_ASF = 0., delta_t_ASF = 0., vs_ASF = 0.;

		if( time_new>=0. and time_new<=t1 ){//C->B
			if( x>0.-eps and x<X1 ){
				delta_m_CB = mB1-mC1;
				delta_c_CB = cB1-cC1;
			}else if( x>X1-eps and x<X2 ){
				delta_m_CB = mB2-mC2;
				delta_c_CB = cB2-cC2;
			}else if( x>X2-eps and x<X3 ){
				delta_m_CB = mB3-mC3;
				delta_c_CB = cB3-cC3;
			}else if( x>X3-eps and x<X4+eps ){
				delta_m_CB = mB4-mC4;
				delta_c_CB = cB4-cC4;
			}else if( x>=X4 and x<=X5 ){
				delta_m_CB = mB5-mC5;
				delta_c_CB = cB5-cC5;
			}
			delta_z_CB = delta_m_CB * (x/1000.) + delta_c_CB;
			delta_t_CB = t1-0.;
			vs_CB = delta_z_CB/delta_t_CB;
			sedimentation_depth = vs_CB * (time_new - 0. );
		}else if( time_new>t1 and time_new<=t2 ){// no sedimentation
			if( x>0.-eps and x<X1 ){
				delta_m_CB = mB1-mC1;
				delta_c_CB = cB1-cC1;
			}else if( x>X1-eps and x<X2 ){
				delta_m_CB = mB2-mC2;
				delta_c_CB = cB2-cC2;
			}else if( x>X2-eps and x<X3 ){
				delta_m_CB = mB3-mC3;
				delta_c_CB = cB3-cC3;
			}else if( x>X3-eps and x<X4+eps ){
				delta_m_CB = mB4-mC4;
				delta_c_CB = cB4-cC4;
			}else if( x>=X4 and x<=X5 ){
				delta_m_CB = mB5-mC5;
				delta_c_CB = cB5-cC5;
			}
			delta_z_CB = delta_m_CB * (x/1000.) + delta_c_CB;
			sedimentation_depth = delta_z_CB;
		}else if( time_new>t2 and time_new<=t3 ){//B->A
			if( x>0.-eps and x<X1 ){
				delta_m_CB = mB1-mC1;
				delta_c_CB = cB1-cC1;
				delta_m_BA = mA1-mB1;
				delta_c_BA = cA1-cB1;
			}else if( x>X1-eps and x<X2 ){
				delta_m_CB = mB2-mC2;
				delta_c_CB = cB2-cC2;
				delta_m_BA = mA2-mB2;
				delta_c_BA = cA2-cB2;
			}else if( x>X2-eps and x<X3 ){
				delta_m_CB = mB3-mC3;
				delta_c_CB = cB3-cC3;
				delta_m_BA = mA3-mB3;
				delta_c_BA = cA3-cB3;
			}else if( x>X3-eps and x<X4+eps ){
				delta_m_CB = mB4-mC4;
				delta_c_CB = cB4-cC4;
				delta_m_BA = mA4-mB4;
				delta_c_BA = cA4-cB4;
			}else if( x>=X4 and x<=X5 ){
				delta_m_CB = mB5-mC5;
				delta_c_CB = cB5-cC5;
				delta_m_BA = mB5-mC5;
				delta_c_BA = cB5-cC5;
			}
			delta_z_CB = delta_m_CB * (x/1000.) + delta_c_CB;
			double sedimentation_depth_CB = delta_z_CB;
			delta_z_BA = delta_m_BA * (x/1000.) + delta_c_BA;
			delta_t_BA = t3-t2;
			vs_BA = delta_z_BA/delta_t_BA;
			double sedimentation_depth_BA = vs_BA * (time_new - t2 );
			sedimentation_depth = sedimentation_depth_CB + sedimentation_depth_BA;
		}else if( time_new>t3 and time_new<=t4 ){// no sedimentation
			if( x>0.-eps and x<X1 ){
				delta_m_CB = mB1-mC1;
				delta_c_CB = cB1-cC1;
				delta_m_BA = mA1-mB1;
				delta_c_BA = cA1-cB1;
			}else if( x>X1-eps and x<X2 ){
				delta_m_CB = mB2-mC2;
				delta_c_CB = cB2-cC2;
				delta_m_BA = mA2-mB2;
				delta_c_BA = cA2-cB2;
			}else if( x>X2-eps and x<X3 ){
				delta_m_CB = mB3-mC3;
				delta_c_CB = cB3-cC3;
				delta_m_BA = mA3-mB3;
				delta_c_BA = cA3-cB3;
			}else if( x>X3-eps and x<X4+eps ){
				delta_m_CB = mB4-mC4;
				delta_c_CB = cB4-cC4;
				delta_m_BA = mA4-mB4;
				delta_c_BA = cA4-cB4;
			}else if( x>=X4 and x<=X5 ){
				delta_m_CB = mB5-mC5;
				delta_c_CB = cB5-cC5;
				delta_m_BA = mB5-mC5;
				delta_c_BA = cB5-cC5;
			}
			delta_z_CB = delta_m_CB * (x/1000.) + delta_c_CB;
			double sedimentation_depth_CB = delta_z_CB;
			delta_z_BA = delta_m_BA * (x/1000.) + delta_c_BA;
			double sedimentation_depth_BA = delta_z_BA;
			sedimentation_depth = sedimentation_depth_CB + sedimentation_depth_BA;
		}else if( time_new>t4 and time_new<=t5 ){//A->SF
			if( x>0.-eps and x<X1 ){
				delta_m_CB = mB1-mC1;
				delta_c_CB = cB1-cC1;
				delta_m_BA = mA1-mB1;
				delta_c_BA = cA1-cB1;
				delta_m_ASF = mSF1-mA1;
				delta_c_ASF = cSF1-cA1;
			}else if( x>X1-eps and x<X2 ){
				delta_m_CB = mB2-mC2;
				delta_c_CB = cB2-cC2;
				delta_m_BA = mA2-mB2;
				delta_c_BA = cA2-cB2;
				delta_m_ASF = mSF2-mA2;
				delta_c_ASF = cSF2-cA2;
			}else if( x>X2-eps and x<X3 ){
				delta_m_CB = mB3-mC3;
				delta_c_CB = cB3-cC3;
				delta_m_BA = mA3-mB3;
				delta_c_BA = cA3-cB3;
				delta_m_ASF = mSF3-mA3;
				delta_c_ASF = cSF3-cA3;
			}else if( x>X3-eps and x<X4+eps ){
				delta_m_CB = mB4-mC4;
				delta_c_CB = cB4-cC4;
				delta_m_BA = mA4-mB4;
				delta_c_BA = cA4-cB4;
				delta_m_ASF = mSF4-mA4;
				delta_c_ASF = cSF4-cA4;
			}else if( x>=X4 and x<=X5 ){
				delta_m_CB = mB5-mC5;
				delta_c_CB = cB5-cC5;
				delta_m_BA = mB5-mC5;
				delta_c_BA = cB5-cC5;
				delta_m_ASF = mSF5-mA5;
				delta_c_ASF = cSF5-cA5;
			}
			delta_z_CB = delta_m_CB * (x/1000.) + delta_c_CB;
			double sedimentation_depth_CB = delta_z_CB;
			delta_z_BA = delta_m_BA * (x/1000.) + delta_c_BA;
			double sedimentation_depth_BA = delta_z_BA;
			delta_z_ASF = delta_m_ASF * (x/1000.) + delta_c_ASF;
			delta_t_ASF = t5-t4;
			vs_ASF = delta_z_ASF/delta_t_ASF;
			double sedimentation_depth_ASF = vs_ASF * (time_new - t4 );
			sedimentation_depth = sedimentation_depth_CB + sedimentation_depth_BA + sedimentation_depth_ASF;
		}

		return sedimentation_depth; /*m*/
	}

	// COMPACTION
	double CompactionFunction( Dune::FieldVector<double,dim> globalpos /*ndim,ndim*/ ) const {

		double x = globalpos[0]*X_c.x_c; /*m*/
		double z = globalpos[dim-1]*X_c.x_c; /*m*/
		double sedimentation_depth /*m*/ = SedimentationDepth(x/*m*/);
		double beta = 1./3000.;
		double zPSFC /*ndim*/ = mesh.z_PSFC(globalpos[0]/*ndim*/);
		zPSFC *= X_c.x_c; /*m*/
		double compaction_factor = std::exp( -beta*(zPSFC-z + sedimentation_depth) );

		return compaction_factor;
	}
	
	// SEDIMENT DEFORMTION/FLOW VELOCITY
	Dune::FieldVector<double,dim>
	SedimentVelocity () const {
//		double time_new = ((*time)+(*dt))*X_c.t_c; /*s*/
		Dune::FieldVector<double,dim> vs( 0. );
		vs[ dim - 1 ] = 0.;
		vs[0] = 0.;

		return vs; /*m/s*/
	}
	
	
	/* GRAVITY VECTOR */
	Dune::FieldVector<double,dim>
	g( ) const {
		Dune::FieldVector<double,dim> gravity( 0. );
		double g = 0.;
		if(gravity_flag) g = gravity_magnitude;
		gravity[ dim - 1 ] = g;
		gravity[0] = 0.;
		return gravity; /*N/kg*/
	}


	/* REFERENCE STATE VALUES */
	double ReferenceSalinity() const {
		return 0.; /*kg/kg*/
	}
	double ReferenceTemperature() const {
		return 273.15/X_c.T_c; /*ndim*/
	}
	double ReferencePressure() const {
		return 1.01e5/X_c.P_c; /*ndim*/
	}

};
