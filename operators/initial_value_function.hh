/** \brief A function for initial values of Pw
 */
template<typename GV, typename Properties, typename PGMap>
class Pw_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Pw_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  Pw_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
                        const Dune::FieldVector<double,dim>& xlocal,
                        double& y) const
  {
//    auto x = element.geometry().global(xlocal);

    y /*ndim*/ = icvalue.evaluate(element,xlocal)[Indices::PVId_Pw] ;
//    std::cout<< "Pw_boundary = " << y << std::endl;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Sg
 */
template<typename GV, typename Properties, typename PGMap>
class Sg_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Sg_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  Sg_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
//    auto x = element.geometry().global(xlocal);

	y = icvalue.evaluate(element,xlocal)[Indices::PVId_Sg] ;
//	std::cout<< "Sg_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};

/** \brief A function for initial values of XCH4
 */
template<typename GV, typename Properties, typename PGMap>
class XCH4_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, XCH4_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  XCH4_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
//    auto x = element.geometry().global(xlocal);

	y = icvalue.evaluate(element,xlocal)[Indices::PVId_XCH4] ;
//	std::cout<< "XCH4_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of YH2O
 */
template<typename GV, typename Properties, typename PGMap>
class YH2O_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, YH2O_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  YH2O_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
//    auto x = element.geometry().global(xlocal);

	y = icvalue.evaluate(element,xlocal)[Indices::PVId_YH2O] ;
//	std::cout<< "YH2O_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Sh
 */
template<typename GV, typename Properties, typename PGMap>
class Sh_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Sh_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  Sh_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
//    auto x = element.geometry().global(xlocal);

    y=icvalue.evaluate(element,xlocal)[Indices::PVId_Sh] ; //initial hydrate saturation
//    std::cout<< "Sh_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of Xc
 */
template<typename GV, typename Properties, typename PGMap>
class Xc_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, Xc_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  Xc_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const
  {
//    auto x = element.geometry().global(xlocal);

    y=icvalue.evaluate(element,xlocal)[Indices::PVId_Xc] ; //initial salt mole fraction
//    std::cout<< "Xc_boundary = " << y << std::endl;

    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of T
 */
template<typename GV, typename Properties, typename PGMap>
class T_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, T_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const PGMap& pg ;
	  const Properties& property;
	  const static int dim = GV::dimension;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;

public:

  //! construct from grid view
  T_Initial ( const GV& gv_, const Properties& property_, const PGMap& pg_ )
  : gv( gv_ ) , property(property_), pg( pg_ ), icvalue(gv_,property_,pg_)
  {}

  //! evaluate extended function on element
  inline void evaluate ( const typename GV::Traits::template Codim<0>::Entity& element,
          	  	  	  	 const Dune::FieldVector<double,dim>& xlocal,
						 double& y) const {

//    auto x = element.geometry().global(xlocal);

    y/*ndim*/ = icvalue.evaluate(element,xlocal)[Indices::PVId_T] ; //initial temperature
//    std::cout<< "T_boundary = " << y << std::endl;

    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
