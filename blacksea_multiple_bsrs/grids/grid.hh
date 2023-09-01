template<typename PTree>
class MeshParameters{
private:
	const PTree& ptree;
	constexpr static double eps = 1.0e-6;
	CharacteristicValues Xc;
	
	double Zmax;
	double Xmax;

public:
	
	//! constructor
	MeshParameters (const PTree& ptree_)
	:ptree(ptree_)
	{
		Zmax=ptree.get("grid.ug.LZ",(double)1000.); /*m*/
		Xmax=ptree.get("grid.ug.LX",(double)1000.); /*m*/
	}

	// (pseudo-1D) 2D -> X and Z
	const static int dimension = 2;
	double Z_length = ptree.get("grid.ug.LZ",(double)1000.)/Xc.x_c; /*ndim*/
	double X_length = ptree.get("grid.ug.LX",(double)1000.)/Xc.x_c; /*ndim*/

	double volumeFactor( double r /*radial component of cell center*/ ) const {

		/* Volume = 2 PI r_center dr dz
		 * */
		double factor = 1.;//
		return factor; /*m^3*/
	}

	double areaFactor( double r_c /*radial component of cell center*/,
					   double r_f /*radial component of face center*/,
					   Dune::FieldVector<double,dimension> normal ) const {

		/* Surface area = 2 PI ( r_center dr i_hat + r_face dz j_hat )
		 * */
		double factor = 1.;
		return abs(factor); /*m^2*/
	}


	bool isBottomBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( globalPos[1] < 0. + eps ){
			return true;
		}
		else
			return false;

	}

	bool isTopBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( globalPos[1] > Zmax/Xc.x_c - eps ){
			return true;
		}
		else
			return false;
	}

};
