
#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "nemo-datastructs.h"
#include "nemo/nemo.h"

// Calculate optical depth to the inner edge of every shell.
void spherical_optical_depth(
    class nemo &,                  ///< pointer to NEMO class
    const nemo_star_properties &,  ///< struct containing stellar properties
    const std::vector<std::vector<double> > &,  ///< vector of state vectors
    std::vector<std::vector<struct nemo_cell_rt_data> > &  ///< vec of cell data
);


#endif  // COMMON_FUNCTIONS_H
