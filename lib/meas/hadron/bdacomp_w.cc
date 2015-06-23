// $Id: bdacomp_w.cc, v 1.0 2011-04-11 chagen, rwschiel $

/*
 *  constructs baryon propagtor. at the source the diquark, constructed
 *  from quark_propagator_1 and quark_propagator_2,
 *  gets contracted with some gamma matrix which can be chosen
 *  from the Clifford algebra via a parameter. the other 4 spin 
 *  indices are left open.
 */

#include "meas/hadron/bdacomp_w.h"

namespace Chroma 
{

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
	        int t0, int bc_spec, int gam_diquark)
  {
    START_CODE ();

    // Length of lattice in decay direction
    int length  = phases.numSubsets ();

    int num_mom = phases.numMom ();

    LatticePropagator q_prop_gamma_1  = quark_propagator_1 * Gamma (gam_diquark);

    multi2d<DComplex> foo (num_mom, length);
    for (int si1 = 0; si1 < Ns; si1++)
      for (int si2 = 0; si2 < Ns; si2++)
        for (int si3 = 0; si3 < Ns; si3++)
          for (int so3 = 0; so3 < Ns; so3++)
            {
              LatticeComplex b_prop = zero;

              for (int s = 0; s < Ns; s++)
              {
                // Contract over color indices with antisym tensors
                b_prop +=
                  colorContract (peekSpin (q_prop_gamma_1, si1, s),
                                 peekSpin (quark_propagator_2, si2, s),
                                 peekSpin (quark_propagator_3, si3, so3));
              }
                foo = phases.sft (b_prop);

                for (int sink_mom_num = 0; sink_mom_num < num_mom; sink_mom_num++) 
                  for (int t = 0; t < length; t++) {
                  //shift source to 0 and take care of the antiperiodic BC
                  //sign flip for the baryon
                    int t_eff = (t - t0 + length) % length;
                    bdaprop [sink_mom_num][t_eff] (si1, si2, si3, so3) = 
                      (bc_spec < 0 && (t_eff + t0) >= length) ? -foo [sink_mom_num][t] :
                       foo [sink_mom_num][t] ;
                  }
            }

    END_CODE ();
  }


  void write_bda (BinaryFileWriter& to,
                  multi2d<BaryonDA>& bdaprop, 
                  const SftMom& phases)
  {

      for(int p = 0; p < phases.numMom (); p++)
      {
        multi1d<DComplex> data (phases.numSubsets () * bdaprop [p][0].handle ().size ());  
        int k = 0;
        for (int t = 0; t < phases.numSubsets (); t++)
          for (int j = 0; j < bdaprop [p][t].Size (); j++)
            data [k++] = bdaprop [p][t][j];

        write (to,data);
      }

   }

}  // end namespace Chroma
