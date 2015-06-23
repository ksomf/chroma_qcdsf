// $Id: bdacomp_w.h, v 1.0 2011-04-11 chagen, rwschiel $

/*
 *  constructs baryon propagtor. at the source the diquark, constructed
 *  from quark_propagator_1 and quark_propagator_2, 
 *  gets contracted with some gamma matrix which can be chosen
 *  from the Clifford algebra via a parameter. the other 4 spin 
 *  indices are left open.
 */

#ifndef __bdacomp_w_h__
#define __bdacomp_w_h__

#include "chromabase.h"
#include "io/qprop_io.h"
#include "util/ft/sftmom.h"

namespace Chroma 
{

  class BaryonDA {
  private:
    multi1d<DComplex> data;
    int size;

  public:

    BaryonDA () {
      size = Ns * Ns * Ns * Ns;
      data.resize (size);
    };

    ~BaryonDA () {};

    DComplex& operator() (int s1, int s2, int s3, int s4) {
      return data [s1 + Ns * (s2 + Ns * (s3 + Ns * s4))];
    };

    multi1d<DComplex>& handle () {
      return data;
    }

    DComplex& operator[] (int i) {
      return data [i];
    }

    int Size () const {return size;}

  };

  //! Construct all components of a baryon propagator
  /*!
   * \ingroup hadron
   *
   * This routine is specific to Wilson fermions!
   *
   * In all baryons the colour components are contracted with the totally
   * antisymmetric 'tensor' eps(a,b,c) = antisym_tensor(a,b,c).
   *
   * \param barprop                  baryon correlation function (in real space) ( Write )
   * \param quark_propagator_1       quark propagator ( Read )
   * \param quark_propagator_2       quark propagator ( Read )
   * \param quark_propagator_3       quark propagator ( Read )
   * \param phases                   object holds list of momenta ( Read )
   * \param t0                       coordinates of source in decay direction ( Read )
   * \param bc_spec                  boundary condition for spectroscopy ( Read )
   */

  void bdacomp (multi2d<BaryonDA>& bdaprop, 
	        const LatticePropagator& quark_propagator_1, 
	        const LatticePropagator& quark_propagator_2,
	        const LatticePropagator& quark_propagator_3,
	        const SftMom& phases,
	        int t0, int bc_spec, int gam_diquark);


  void write_bda (BinaryFileWriter& to,
                  multi2d<BaryonDA>& bdaprop, 
                  const SftMom& phases);

}  // end namespace Chroma

#endif
