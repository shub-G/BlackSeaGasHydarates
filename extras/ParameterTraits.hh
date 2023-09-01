/*
 * ParameterTraits.hh
 *
 *  Created on: Nov 16, 2016
 *      Author: shubhangi
 */

#ifndef EXTRAS_PARAMETERTRAITS_HH_
#define EXTRAS_PARAMETERTRAITS_HH_

template<typename GV, typename RF>
struct ParameterTraits
{
  //! \brief the grid view
  typedef GV GridViewType;

  //! \brief Enum for domain dimension
  enum {
    //! \brief dimension of the domain
    dimDomain = GV::dimension
  };

  //! \brief Export type for domain field
  typedef typename GV::Grid::ctype DomainFieldType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

  //! \brief domain type
  typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

  //! \brief Export type for range field
  typedef RF RangeFieldType;

  //! \brief range type
  typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

  //! \brief permeability tensor type
  typedef Dune::FieldMatrix<RangeFieldType,dimDomain,dimDomain> PermTensorType;

  //! grid types
  typedef typename GV::Traits::template Codim<0>::Entity ElementType;
  typedef typename GV::Intersection IntersectionType;
};



#endif /* EXTRAS_PARAMETERTRAITS_HH_ */
