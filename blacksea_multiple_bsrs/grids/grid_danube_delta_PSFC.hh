/*
 * grid_danube_delta_PSFC.hh
 *
 *  Created on: Apr 12, 2022
 *      Author: sgupta
 */

#ifndef BLACKSEA_MULTIPLE_BSRS_GRIDS_GRID_DANUBE_DELTA_PSFC_HH_
#define BLACKSEA_MULTIPLE_BSRS_GRIDS_GRID_DANUBE_DELTA_PSFC_HH_

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
		Zmax=ptree.get("grid.ug.LZ",(double)800.); /*m*/
		Zmax *= 1./Xc.x_c; /*ndim*/
		Xmax=ptree.get("grid.ug.LX",(double)20000.); /*m*/
		Xmax *= 1./Xc.x_c; /*ndim*/
	}

	/*
	 * 2D -> X and Z
	 */
	const static int dimension 		= 2  	;
	constexpr static double origin 		= 0.0	;
#ifdef RADIAL
	constexpr static double well_radius 	= 1.0/Xc.x_c		;
#endif

	double volumeFactor( double r /*radial component of cell center*/ ) const {

		/* Volume = 2 PI r_center dr dz
		 * */
		double factor = 1.;//
#ifdef RADIAL
		r += well_radius;
		factor = 2.*PI*r;
#endif
		return factor;
	}

	double areaFactor( double r_c /*radial component of cell center*/,
					   double r_f /*radial component of face center*/,
					   Dune::FieldVector<double,dimension> normal ) const {

		/* Surface area = 2 PI ( r_center dr i_hat + r_face dz j_hat )
		 * */
		double factor = 1.;
#ifdef RADIAL
		r_c += well_radius;
		r_f += well_radius;
		factor = 2.*PI*(r_c*normal[0]+r_f*normal[1]);
#endif
		return std::abs(factor);
	}


	bool isLeftBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( globalPos[0] < origin + eps )
			return true;
		else return false;
	}

	bool isRightBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( globalPos[0] > Xmax - eps )
			return true;
		else return false;
	}

	bool isBottomBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( globalPos[1] < origin + eps ){
			return true;
		}
		else return false;
	}

	bool isTopBoundary( Dune::FieldVector< double, dimension > globalPos /*ndim*/ ) const{
		if( !isLeftBoundary(globalPos) and !isRightBoundary(globalPos) and !isBottomBoundary(globalPos) ){
			return true;
		}
		else return false;
	}

	double z_PSFC( double xpos /* global x-location, ndim */ )const {
		int num_of_points = 26;
		std::vector<double> x(num_of_points,0.);
		std::vector<double> y(num_of_points,0.);
		double x_c = Xc.x_c;
		x[0]  = 0*1.e3/x_c		, y[0]  = 30.0/x_c	;
		x[1]  = 1.*1.e3/x_c		, y[1]  = 45.0/x_c	;
		x[2]  = 2.*1.e3/x_c		, y[2]  = 70.0/x_c	;
		x[3]  = 3.*1.e3/x_c		, y[3]  = 80.0/x_c	;
		x[4]  = 4.*1.e3/x_c		, y[4]  = 95.0/x_c	;
		x[5]  = 5.*1.e3/x_c		, y[5]  = 120.0/x_c	;
		x[6]  = 6.*1.e3/x_c		, y[6]  = 130.0/x_c	;
		x[7]  = 7.*1.e3/x_c		, y[7]  = 180.0/x_c	;
		x[8]  = 8.*1.e3/x_c		, y[8]  = 245.0/x_c	;
		x[9]  = 8.5*1.e3/x_c	, y[9]  = 260.0/x_c	;
		x[10] = 9.*1.e3/x_c		, y[10] = 235.0/x_c	;
		x[11] = 9.5*1.e3/x_c	, y[11] = 185.0/x_c	;
		x[12] = 10.*1.e3/x_c	, y[12] = 130.0/x_c	;
		x[13] = 10.5*1.e3/x_c	, y[13] = 100.0/x_c	;
		x[14] = 11.*1.e3/x_c	, y[14] = 85.0/x_c	;
		x[15] = 11.5*1.e3/x_c	, y[15] = 75.0/x_c	;
		x[16] = 12.*1.e3/x_c	, y[16] = 105.0/x_c	;
		x[17] = 12.5*1.e3/x_c	, y[17] = 185.0/x_c	;
		x[18] = 13.*1.e3/x_c	, y[18] = 210.0/x_c	;
		x[19] = 14.*1.e3/x_c	, y[19] = 145.0/x_c	;
		x[20] = 15.*1.e3/x_c	, y[20] = 110.0/x_c	;
		x[21] = 16.*1.e3/x_c	, y[21] = 105.0/x_c	;
		x[22] = 17.*1.e3/x_c	, y[22] = 75.0/x_c	;
		x[23] = 18.*1.e3/x_c	, y[23] = 55.0/x_c	;
		x[24] = 19.*1.e3/x_c	, y[24] = 10.0/x_c	;
		x[25] = 20.*1.e3/x_c	, y[25] = 0.0/x_c	;


		double zPSFC = 0.;
		for( int i=0; i<num_of_points-1; i++ ){
			if( xpos > x[i] and xpos < x[i+1]+1.e-6){
				zPSFC = y[i] + (y[i+1]-y[i])/(x[i+1]-x[i])*(xpos-x[i]);
				i = num_of_points;
			}else{
				zPSFC = 0.;
			}
		}

//		std::cout << xpos << '\t' << zPSFC << std::endl;
		return Zmax+zPSFC; /*ndim*/
	}
};

#endif /* BLACKSEA_MULTIPLE_BSRS_GRIDS_GRID_DANUBE_DELTA_PSFC_HH_ */
