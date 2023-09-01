namespace Dune {
  namespace PDELab {


  template<class GFS, class Vector>
  class Evaluation{
  public:

	  double dt;

	  typedef typename GFS::Traits::GridViewType GV;
	  typedef LocalFunctionSpace<GFS> LFS;
	  typedef Dune::PDELab::LFSIndexCache<LFS> LFSCache;
	  typedef typename Vector::template LocalView<LFSCache> VectorView;
	  typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType RangeType;
      typedef typename LFS::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType JacobianType;
      typedef ParameterTraits<GV,double> Traits;


	  Evaluation(GFS gfs_,Vector *xglobal_):
		  gfs(gfs_),xglobal(xglobal_),xglobal_old(xglobal_),xglobal_fp(xglobal_),dtxglobal(xglobal_),dt(1.){
	  }

	  Evaluation(GFS gfs_,Vector *xglobal_,Vector *xglobal_old_):
		  gfs(gfs_),xglobal(xglobal_),xglobal_old(xglobal_old_),xglobal_fp(xglobal_),dtxglobal(xglobal_){
		  dt=1e-5;
	  }

	  Evaluation(GFS gfs_,Vector *xglobal_,Vector *xglobal_old_, double dt_):
		  gfs(gfs_),xglobal(xglobal_),xglobal_old(xglobal_old_),xglobal_fp(xglobal_),dtxglobal(xglobal_),dt(dt_){
	  }

	  Evaluation(GFS gfs_,Vector *xglobal_,Vector *xglobal_old_,Vector *xglobal_fp_, double dt_):
		  gfs(gfs_),xglobal(xglobal_),xglobal_old(xglobal_old_),xglobal_fp(xglobal_fp_),dtxglobal(xglobal_),dt(dt_){
	  }


	  Evaluation(GFS gfs_,Vector *xglobal_,Vector *xglobal_old_,Vector *dtxglobal_,Vector *xglobal_fp_, double dt_):
		  gfs(gfs_),xglobal(xglobal_),xglobal_old(xglobal_old_),xglobal_fp(xglobal_fp_),dtxglobal(dtxglobal_),dt(dt_){
	  }

	  void update(Vector *xglobal_){
		  dtxglobal = *xglobal_old;
		  xglobal_old = xglobal;
		  xglobal=xglobal_;
	  }

	  void update(Vector *xglobal_,Vector *xglobal_old_){
		  xglobal=xglobal_;
		  xglobal_old=xglobal_old_;
	  }

	  void update(Vector *xglobal_,Vector *xglobal_old_, double dt_){
		  xglobal=xglobal_;
		  xglobal_old=xglobal_old_;
		  dt=dt_;
	  }

	  void update(Vector *xglobal_,Vector *xglobal_old_,Vector *xglobal_fp_, double dt_){
		  xglobal=xglobal_;
		  xglobal_old=xglobal_old_;
		  xglobal_fp=xglobal_fp_;
		  dt=dt_;
	  }

	  void update(Vector *xglobal_,Vector *xglobal_old_,Vector *xglobal_fp_,Vector *dtxglobal_, double dt_){
		  xglobal=xglobal_;
		  xglobal_old=xglobal_old_;
		  xglobal_fp=xglobal_fp_;
		  dtxglobal=dtxglobal_;
		  dt=dt_;
	  }

	  void update(double dt_){
		  dt=dt_;
	  }

	  void updateFP(Vector *xglobal_fp_){
		  xglobal_fp=xglobal_fp_;
	  }


	  inline void evalFunction(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = 0;
		  for (unsigned int i=0; i<xlocal.size(); i++){
			  (*y) += xlocal[i]*phi[i];
		  }

	  }

	  inline void evalFunctionMass(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = xlocal[0]*phi[0];

	  }


	  inline void evalGradient(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,Dune::FieldVector<double,GV::dimension>  *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space
		  const int dim = GV::dimension;

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  //
		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);

		  // get local Jacobians/gradients of the shape functions
		  std::vector<JacobianType> J(lfs.size());
		  lfs.finiteElement().localBasis().evaluateJacobian(x,J);

		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  jac = e.geometry().jacobianInverseTransposed(x);

		  // Initialize vectors
		  std::vector<Dune::FieldVector<RangeType,GV::dimension> > gradphi(lfs.size());
		  for (unsigned int  i=0; i<lfs.size(); i++)
			  jac.mv(J[i][0],gradphi[i]);

		  (*y)=0.;
		  Dune::FieldVector<RangeType,GV::dimension> gradu(0.0);
		  for (unsigned int i=0; i<lfs.size(); i++)
			  gradu.axpy(xlocal[i],gradphi[i]);
		  (*y)=gradu;

//		  //typename Traits::RangeType gradphi;
//		  (*y) = 0.;
//		  Dune::FieldVector<RangeType,GV::dimension> gradphi(0.0);
//
//          std::vector<Dune::FieldVector<RangeType,GV::dimension> > tgradphi_s(lfs.size());
//          for (unsigned int i=0; i<lfs.size(); i++) JgeoIT.mv(J[i][0],tgradphi_s[i]);
//
////          std::vector<Dune::FieldVector<RangeType,GV::dimension> > gradphi(lfs.size());
////          for (int i=0; i<lfs.size(); i++)
////        	  JgeoIT.mv(J[i][0],gradphi[i]);
////
////          // compute gradient of u
////          Dune::FieldVector<RF,GV::dimension> gradu(0.0);
////          for (int i=0; i<lfs.size(); i++)
////            gradu.axpy(xl(lfs,i),gradphi[i]);
////
////          (*y)=gradu;
//
//          Dune::FieldVector<RangeType,GV::dimension> gradu(0.0);
//          for (unsigned int i=0; i<lfs.size(); i++){
//        	  gradphi = 0;
//        	  gradphi = tgradphi_s[i];
//        	  gradphi*=xlocal[i];
//        	  gradu+=gradphi;
//          }
//
//          (*y)=gradu;
//
////		  for(unsigned i = 0; i < lfs.size(); ++i) {
////			  // compute global gradient of shape function i
////			  gradphi = 0;
////			  JgeoIT.umv(J[i][0], gradphi);
////			  gradphi*=xlocal[i];
////			  (*y)+=gradphi;
////			  // sum up global gradients, weighting them with the appropriate coeff
////			  //grad.axpy(xlocal[i], gradphi);
////		  }
//		  //y=grad;
	  }


//	  inline void evalGradientWeighted(const typename Traits::ElementType& e,
//			  const typename Traits::DomainType& x,Dune::FieldVector<double,GV::dimension>  *y){
//
//		  const int dim = GV::dimension;
//
//		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
//		  // get and bind local functions space
//		  LFS lfs(gfs);
//		  LFSCache lfs_cache(lfs);
//		  VectorView x_view((*xglobal));
//
//		  lfs.bind(e);
//		  lfs_cache.update();
//		  x_view.bind(lfs_cache);
//
//		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
//		  x_view.read(xlocal);
//
//		  // get Jacobian of geometry
//		  const typename Traits::ElementType::Geometry::Jacobian&
//		  JgeoIT = e.geometry().jacobianInverseTransposed(x);
//
//          typename Traits::PermTensorType A;
//          A = parameter.A(e.entity(),x);
//
//		  std::vector<RangeType> phi(lfs.size(),0.);
//		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
//		  double u = 0;
//		  for (unsigned int i=0; i<xlocal.size(); i++){
//			  u += xlocal[i]*phi[i];
//		  }
//          double Hu=parameter.H(u);
//
//		  // get local Jacobians/gradients of the shape functions
//		  std::vector<JacobianType> J(lfs.size());
//		  lfs.finiteElement().localBasis().evaluateJacobian(x,J);
//		  //typename Traits::RangeType gradphi;
//		  (*y) = 0.;
//		  Dune::FieldVector<RangeType,GV::dimension> gradphi(0.0);
//
//          // compute gradient of u
//          Dune::FieldVector<double,GV::dimension> gradu(0.0);
//
//		  for(unsigned i = 0; i < lfs.size(); ++i) {
//			  // compute global gradient of shape function i
//			  gradphi = 0;
//			  JgeoIT.umv(J[i][0], gradphi);
//			  gradphi*=xlocal[i];
//			  gradu+=gradphi;
//			  // sum up global gradients, weighting them with the appropriate coeff
//			  //grad.axpy(xlocal[i], gradphi);
//		  }
//
//          A.umv(gradu,(*y));
//		  //y=grad;
//	  }

	  inline void evalFunctionOld(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = 0;
		  for (unsigned int i=0; i<xlocal.size(); i++){
			  (*y) += xlocal[i]*phi[i];
		  }

	  }

	  inline void evalFunctionOldMass(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = xlocal[0]*phi[0];

	  }

	  inline void evalGradientOld(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,Dune::FieldVector<double,GV::dimension>  *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space
		  const int dim = GV::dimension;

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  JgeoIT = e.geometry().jacobianInverseTransposed(x);

		  // get local Jacobians/gradients of the shape functions
		  std::vector<JacobianType> J(lfs.size());
		  lfs.finiteElement().localBasis().evaluateJacobian(x,J);
		  //typename Traits::RangeType gradphi;
		  (*y) = 0.;
		  Dune::FieldVector<RangeType,GV::dimension> gradphi(0.0);

          std::vector< Dune::FieldVector<RangeType,GV::dimension> > tgradphi_s(lfs.size());
          for (unsigned int i=0; i<lfs.size(); i++) JgeoIT.mv(J[i][0],tgradphi_s[i]);

          Dune::FieldVector<RangeType,GV::dimension> gradu(0.0);
          for (unsigned int i=0; i<lfs.size(); i++){
        	  gradphi = 0;
        	  gradphi = tgradphi_s[i];
        	  gradphi*=xlocal[i];
        	  for(int vec_i=0;vec_i<dim;vec_i++){
            	  gradu[vec_i]+=gradphi[vec_i];
        	  }

          }

          (*y)=gradu;

//		  for(unsigned i = 0; i < lfs.size(); ++i) {
//			  // compute global gradient of shape function i
//			  gradphi = 0;
//			  JgeoIT.umv(J[i][0], gradphi);
//			  gradphi*=xlocal[i];
//			  (*y)+=gradphi;
//			  // sum up global gradients, weighting them with the appropriate coeff
//			  //grad.axpy(xlocal[i], gradphi);
//		  }
		  //y=grad;
	  }

	  inline void evaldtFunction(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));
		  VectorView x_view_old((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);
		  x_view_old.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  std::vector<typename Vector::ElementType> xlocal_old(lfs.size());
		  x_view.read(xlocal);
		  x_view_old.read(xlocal_old);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = 0;
		  for (unsigned int i=0; i<xlocal.size(); i++){
			  (*y) += (xlocal[i]/dt-xlocal_old[i]/dt)*phi[i];
		  }

	  }


	  inline void evaldtFunctionOld(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*dtxglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = 0;
		  for (unsigned int i=0; i<xlocal.size(); i++){
			  (*y) += xlocal[i]*phi[i];
		  }

	  }


	  inline void evaldtFunctionMass(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));
		  VectorView x_view_old((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);
		  x_view_old.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  std::vector<typename Vector::ElementType> xlocal_old(lfs.size());
		  x_view.read(xlocal);
		  x_view_old.read(xlocal_old);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = (xlocal[0]/dt-xlocal_old[0]/dt)*phi[0];

	  }


	  inline void evaldtFunctionOldMass(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*dtxglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = xlocal[0]*phi[0];

	  }

	  inline void evaldtGradient(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,Dune::FieldVector<double,GV::dimension>  *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space
		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));
		  VectorView x_view_old((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);
		  x_view_old.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  std::vector<typename Vector::ElementType> xlocal_old(lfs.size());
		  x_view.read(xlocal);
		  x_view_old.read(xlocal_old);

		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  JgeoIT = e.geometry().jacobianInverseTransposed(x);

		  // get local Jacobians/gradients of the shape functions
		  std::vector<JacobianType> J(lfs.size());
		  lfs.finiteElement().localBasis().evaluateJacobian(x,J);
		  //typename Traits::RangeType gradphi;
		  (*y) = 0.;
		  Dune::FieldVector<RangeType,GV::dimension> gradphi(0.0);

		  for(unsigned i = 0; i < lfs.size(); ++i) {
			  // compute global gradient of shape function i
			  gradphi = 0;
			  JgeoIT.umv(J[i][0], gradphi);
			  gradphi*=(xlocal[i]-xlocal_old[i])/dt;
			  //(*y)+=gradphi;
        	  for(int vec_i=0;vec_i<GV::dimension;vec_i++){
        		  (*y)[vec_i]+=gradphi[vec_i];
        	  }
			  // sum up global gradients, weighting them with the appropriate coeff
			  //grad.axpy(xlocal[i], gradphi);
		  }
		  //y=grad;
	  }


	  inline void evalFunctionFP(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_fp));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = 0;
		  for (unsigned int i=0; i<xlocal.size(); i++){
			  (*y) += xlocal[i]*phi[i];
		  }

	  }

	  inline void evalFunctionFPMass(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,
			  double *y) const{

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_fp));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  std::vector<RangeType> phi(lfs.size(),0.);
		  lfs.finiteElement().localBasis().evaluateFunction(x,phi);
		  (*y) = xlocal[0]*phi[0];

	  }


	  inline void evalGradientFP(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,Dune::FieldVector<double,GV::dimension>  *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space

		  const int dim = GV::dimension;

		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_fp));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  JgeoIT = e.geometry().jacobianInverseTransposed(x);

		  // get local Jacobians/gradients of the shape functions
		  std::vector<JacobianType> J(lfs.size());
		  lfs.finiteElement().localBasis().evaluateJacobian(x,J);
		  //typename Traits::RangeType gradphi;
		  (*y) = 0.;
		  Dune::FieldVector<RangeType,GV::dimension> gradphi(0.0);

          std::vector<Dune::FieldVector<RangeType,GV::dimension> > tgradphi_s(lfs.size());
          for (unsigned int i=0; i<lfs.size(); i++) JgeoIT.mv(J[i][0],tgradphi_s[i]);

//          std::vector<Dune::FieldVector<RangeType,GV::dimension> > gradphi(lfs.size());
//          for (int i=0; i<lfs.size(); i++)
//        	  JgeoIT.mv(J[i][0],gradphi[i]);
//
//          // compute gradient of u
//          Dune::FieldVector<RF,GV::dimension> gradu(0.0);
//          for (int i=0; i<lfs.size(); i++)
//            gradu.axpy(xlocal(lfs,i),gradphi[i]);
//
//          (*y)=gradu;

          Dune::FieldVector<RangeType,GV::dimension> gradu(0.0);
          for (unsigned int i=0; i<lfs.size(); i++){
        	  gradphi = 0;
        	  gradphi = tgradphi_s[i];
        	  gradphi*=xlocal[i];
        	  //gradu+=gradphi;
        	  for(int vec_i=0;vec_i<dim;vec_i++){
        		  gradu[vec_i]+=gradphi[vec_i];
        	  }
          }

          (*y)=gradu;

//		  for(unsigned i = 0; i < lfs.size(); ++i) {
//			  // compute global gradient of shape function i
//			  gradphi = 0;
//			  JgeoIT.umv(J[i][0], gradphi);
//			  gradphi*=xlocal[i];
//			  (*y)+=gradphi;
//			  // sum up global gradients, weighting them with the appropriate coeff
//			  //grad.axpy(xlocal[i], gradphi);
//		  }
		  //y=grad;
	  }



	  inline void evalLaplace(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,double *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space
		  const int dim = GV::dimension;
		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

//		  std::vector< Dune::FieldMatrix<RangeType,GV::dimension,GV::dimension> > HesseMatrices(lfs.size());
//		  HesseMatrices[3][0][0] = 1.;HesseMatrices[3][0][1] = 0.;HesseMatrices[3][1][0] = 0.;HesseMatrices[3][1][1] = 0.;
//		  HesseMatrices[4][0][0] = 0.;HesseMatrices[4][0][1] = 1.;HesseMatrices[4][1][0] = 1.;HesseMatrices[4][1][1] = 0.;
//		  HesseMatrices[5][0][0] = 0.;HesseMatrices[5][0][1] = 0.;HesseMatrices[5][1][0] = 0.;HesseMatrices[5][1][1] = 1.;
		  Dune::FieldMatrix<RangeType,GV::dimension,GV::dimension> HesseMatrix(0.);
		  HesseMatrix[0][0] = xlocal[3];HesseMatrix[0][1] = xlocal[4];HesseMatrix[1][0] = xlocal[4];HesseMatrix[1][1] = xlocal[5];
//		  for(int i=3;i<6;i++){
//			  HesseMatrices[i]*=xlocal[i];
////			  std::cout<<" HesseMatrices = "<<HesseMatrices[i]<<std::endl;
//			  HesseMatrix+=HesseMatrices[i];
//		  }


		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  JgeoIT = e.geometry().jacobianInverseTransposed(x);
		  const typename Traits::ElementType::Geometry::JacobianTransposed&
		  JgeoT = e.geometry().jacobianTransposed(x);

		  //std::cout<<"1. HesseMatrix = "<<HesseMatrix<<std::endl;
		  HesseMatrix=HesseMatrix.leftmultiply(JgeoIT);
		  //std::cout<<"2. HesseMatrix = "<<HesseMatrix<<std::endl;
		  HesseMatrix=HesseMatrix.leftmultiply(JgeoIT);
		  //std::cout<<"3. HesseMatrix = "<<HesseMatrix<<std::endl;

		  (*y) = HesseMatrix[0][0] + HesseMatrix[1][1];

		  //y=grad;
	  }


	  inline void evalLaplaceOld(const typename Traits::ElementType& e,
			  const typename Traits::DomainType& x,double *y){

		  //std::cout<<"Werte Gradient "<<name<<" aus."<<std::endl;
		  // get and bind local functions space
		  const int dim = GV::dimension;
		  LFS lfs(gfs);
		  LFSCache lfs_cache(lfs);
		  VectorView x_view((*xglobal_old));

		  lfs.bind(e);
		  lfs_cache.update();
		  x_view.bind(lfs_cache);

		  std::vector<typename Vector::ElementType> xlocal(lfs.size());
		  x_view.read(xlocal);

//		  std::vector< Dune::FieldMatrix<RangeType,GV::dimension,GV::dimension> > HesseMatrices(lfs.size());
//		  HesseMatrices[3][0][0] = 1.;HesseMatrices[3][0][1] = 0.;HesseMatrices[3][1][0] = 0.;HesseMatrices[3][1][1] = 0.;
//		  HesseMatrices[4][0][0] = 0.;HesseMatrices[4][0][1] = 1.;HesseMatrices[4][1][0] = 1.;HesseMatrices[4][1][1] = 0.;
//		  HesseMatrices[5][0][0] = 0.;HesseMatrices[5][0][1] = 0.;HesseMatrices[5][1][0] = 0.;HesseMatrices[5][1][1] = 1.;
		  Dune::FieldMatrix<RangeType,GV::dimension,GV::dimension> HesseMatrix(0.);
		  HesseMatrix[0][0] = xlocal[3];HesseMatrix[0][1] = xlocal[4];HesseMatrix[1][0] = xlocal[4];HesseMatrix[1][1] = xlocal[5];
//		  for(int i=3;i<6;i++){
//			  HesseMatrices[i]*=xlocal[i];
////			  std::cout<<" HesseMatrices = "<<HesseMatrices[i]<<std::endl;
//			  HesseMatrix+=HesseMatrices[i];
//		  }


		  // get Jacobian of geometry
		  const typename Traits::ElementType::Geometry::JacobianInverseTransposed&
		  JgeoIT = e.geometry().jacobianInverseTransposed(x);
		  const typename Traits::ElementType::Geometry::JacobianTransposed&
		  JgeoT = e.geometry().jacobianTransposed(x);

		  //std::cout<<"1. HesseMatrix = "<<HesseMatrix<<std::endl;
		  HesseMatrix=HesseMatrix.leftmultiply(JgeoIT);
		  //std::cout<<"2. HesseMatrix = "<<HesseMatrix<<std::endl;
		  HesseMatrix=HesseMatrix.leftmultiply(JgeoIT);
		  //std::cout<<"3. HesseMatrix = "<<HesseMatrix<<std::endl;

		  (*y) = HesseMatrix[0][0] + HesseMatrix[1][1];

		  //y=grad;
	  }


  private:
	  GFS gfs;
	  Vector *xglobal;
	  Vector *xglobal_fp;
	  Vector *xglobal_old;
	  Vector *dtxglobal;


  };



}
}
