
#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H

#include "nemo-datastructs.h"
#include "nemo/nemo.h"

// Calculate optical depth to the inner edge of every shell.
void spherical_optical_depth(
    class nemo &,          ///< pointer to NEMO class
    const rad_src_info &,  ///< struct containing stellar properties
    const std::vector<std::vector<double> > &,  ///< vector of state vectors
    const std::vector<double> &,                ///< vector of cell diameters
    std::vector<std::vector<double> >
        &  ///< [output] vector of optical depths to inner edge of shell
);


#endif  // COMMON_FUNCTIONS_H
