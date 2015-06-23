// -*- C++ -*-
/*! \file
 *  \brief Calculates the RMS of a wavefunction
 */

#ifndef __rms_w_h__
#define __rms_w_h__

#include "chromabase.h"

namespace Chroma 
{


  void rms(const LatticePropagator& propagator, const int j_decay, 
	   multi1d<int>& srcloc, const bool psi_wavefn,const bool psi_psi_wavefn,
	   XMLWriter& xml, const string& xml_group);

}  


#endif
