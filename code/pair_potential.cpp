/**
 Pair Potentials
 \author Nathan A. Mahynski
 **/

#include "pair_potential.h"

namespace pair_potential {

/*!
  Set the cutoff values for the potential
  \param [in] rcut Cutoff radius
  \param [in] Ushift Energy shift to apply at r \ge r_{cut}
*/
void Potential::set_cutoff (const double nrcut, const double nUshift) {
  rcut = rcut;
  Ushift = Ushift;
}

/*!
  \param [in] \*coeffs Should be passed as [\epsilon, \sigma, \Delta]
*/
void slj::set_coeff (const vector <double> coeffs) {
  if (coeffs.size() != 3) {
    //flag_error ("Failed to set coefficients for SLJ", __FILE__, __LINE__);
    exit(-1);
  }
  epsilon_ = coeffs[0];
  sigma_ = coeffs[1];
  delta_ = coeffs[2];
}
	
/*!
  \param [in] \*a1 Pointer to atom1
  \param [in] \*a2 Pointer to atom2
  \param [in] r2 Minimum image distance squared between atoms
*/
inline double slj::energy (const double r2) {
  double r = sqrt(r2), x = r - delta_;
  if (x < rcut) {
    double a = sigma_/x, a2 = a*a, a6 = a2*a2*a2;
    return 4.0*epsilon_*(a6*a6-a6)+Ushift;
  } else {
    return 0.0;
  }
}
	
/*!
  The vector xyz MUST be pointing from a1 to a2 for the vectors to be correct.  This is automatically handled in min_image_dist2().
  \param [in] r2 Minimum image distance squared between atoms
  \param [in] \*xyz The minimum image vector from a1 to a2
*/
inline vector<double> slj::force (const double r2, const double *xyz) {
  double r = sqrt(r2), x = r - delta_;
  vector <double> force_vec(3,0.0);
  if (x < rcut) {
    double b = 1.0/x, a = sigma_*b, a2 = a*a, a6 = a2*a2*a2, val, factor;
    factor = 24.0*epsilon_*a6*(2.0*a6-1.0)*b/r;
    for (int i = 0; i < 3; ++i) {
      //val = xyz[i]*factor;
      force_vec[i] = xyz[i]*factor;
      //a1->force[i] -= val;
      //a2->force[i] += val;
    }
  } 
  return force_vec;
}

}
