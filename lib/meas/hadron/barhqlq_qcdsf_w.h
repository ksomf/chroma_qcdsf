// -*- C++ -*-
// $Id: barhqlq_w.h,v 3.2 2007-11-30 06:38:20 kostas Exp $
/*! \file
 *  \brief Heavy-light baryon 2-pt functions
 */

#ifndef __barhqlq_qcdsf_w_h__
#define __barhqlq_qcdsf_w_h__

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/barhqlq_w.h"

namespace Chroma 
{


  struct Baryons_corr_QCDSF_t
  {
    //    multi1d<int>     insert_mom;
    multi1d<Complex>       correlator;
  };
  
  struct Baryons_N_QCDSF_t
  {
    multi1d<Baryons_corr_QCDSF_t >     momentum;
  };

  struct BaryonsQCDSF_t
  {
    multi1d<Baryons_N_QCDSF_t > barnum;
  };



  void write(BinaryWriter& bin, const Baryons_corr_QCDSF_t& mes);
  void write(BinaryWriter& bin, const Baryons_N_QCDSF_t& mes);
  void write(BinaryWriter& bin, const BaryonsQCDSF_t& baryons);



  void barhqlq_qcdsf_xml(const LatticePropagator& propagator_1, 
			 const LatticePropagator& propagator_2, 
			 const LatticePropagator& propagator_3,
			 const bool & haveThird,
			 const SftMom& phases,
			 int t0, int bc_spec, bool time_rev, bool fwdbwd_average,
			 XMLWriter& xml,
			 const string& xml_group);


  void barhqlq_qcdsf_lime(const LatticePropagator& propagator_1, 
			  const LatticePropagator& propagator_2, 
			  const LatticePropagator& propagator_3,
			  const bool & haveThird,
			  const SftMom& phases,
			  int t0, int bc_spec, bool time_rev, bool fwdbwd_average,
			  BaryonsQCDSF_t& bar,
			  BaryonsQCDSF_t& bar_trev);

  void barhqlq_qcdsf(const LatticePropagator& propagator_1,
		     const LatticePropagator& propagator_2,
		     const LatticePropagator& propagator_3,
		     const bool & haveThird,
		     const SftMom& phases,
		     multi3d<DComplex>& barprop);

  
}


#endif
