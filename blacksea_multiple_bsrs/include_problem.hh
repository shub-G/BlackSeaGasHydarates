#ifndef BLACKSEA_MULTIPLE_BSRS_INCLUDE_PROBLEM_HH_
#define BLACKSEA_MULTIPLE_BSRS_INCLUDE_PROBLEM_HH_

#include "../extras/ParameterTraits.hh"
#include "../extras/Evaluation.hh"

#include "../operators/indices.hh"

/*******************************************/
// PROBLEM SPECIFICATION
#include"characteristic_values.hh"
#include"grids/grid_danube_delta_PSFC.hh"
#include"parameters.hh"
#include"properties/include_properties.hh"
#include"initial_conditions.hh"
#include"boundary_conditions.hh"
/*******************************************/
#ifdef PLOT_VELOCITIES
#include "../operators/phase_velocity.hh"
#include "../operators/operator_l2projection.hh"
#endif
#include "../operators/post_process.hh"
#include "../operators/initial_value_function.hh"
#include "../operators/operator_local_new.hh"
#include "../operators/operator_time_new.hh"
#include "driver.hh"

#endif /* BLACKSEA_MULTIPLE_BSRS_INCLUDE_PROBLEM_HH_ */
