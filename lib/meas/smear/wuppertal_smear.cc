// $Id: wuppertal_smear.cc, v 1.0 2011-04-08 ehmann, rwschiel $
/*! \file
 *  \brief Wuppertal smearing of color vector
 */

#include "chromabase.h"
#include "meas/smear/wuppertal_smear.h"
#include "actions/boson/operator/klein_gord.h"

namespace Chroma 
{

  //! Do a covariant Wuppertal smearing of a lattice field
  /*!
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param kappa    hopping parameter ( Read )
   *  \param ItrJacobi  number of iterations ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  template<typename T>
  void wuppertalSmear (const multi1d<LatticeColorMatrix>& u, 
		       T& chi, 
		       const Real& kappa, int ItrJacobi, int j_decay)
  {
    T psi;
    Real norm;

     if (j_decay < Nd)
      norm = Real (1) + Real (2 * Nd - 2) * kappa;
     else
      norm = Real (1) + Real (2 * Nd) * kappa;
    
    for (int n = 0; n < ItrJacobi; n++)
    {
      psi = chi;

     for (int mu = 0; mu < Nd; mu++)
      if (mu != j_decay)
      {
	chi += kappa * (u [mu] * shift (psi, FORWARD, mu) + shift (adj (u [mu]) * psi, BACKWARD, mu));
      }
      chi = chi / norm; 
    } 
  }

  //! Do a covariant Wuppertal smearing of a lattice color vector field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      color vector field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void wuppertalSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeColorVector& chi, 
		       const Real& width, int ItrGaus, int j_decay)
  {
    wuppertalSmear<LatticeColorVector> (u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Wuppertal smearing of a lattice fermion field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      fermion field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void wuppertalSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeFermion& chi, 
		       const Real& width, int ItrGaus, int j_decay)
  {
    wuppertalSmear<LatticeFermion> (u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Wuppertal smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      propagator field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void wuppertalSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeStaggeredPropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay)
  {
    wuppertalSmear<LatticeStaggeredPropagator> (u, chi, width, ItrGaus, j_decay);
  }


  //! Do a covariant Wuppertal smearing of a lattice propagator field
  /*! This is a wrapper over the template definition
   *
   * \ingroup smear
   *
   * Arguments:
   *
   *  \param u        gauge field ( Read )
   *  \param chi      propagator field ( Modify )
   *  \param width    width of "shell" wave function ( Read )
   *  \param ItrGaus  number of iterations to approximate Gaussian ( Read )
   *  \param j_decay  direction of decay ( Read )
   */

  void wuppertalSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticePropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay)
  {
    wuppertalSmear<LatticePropagator> (u, chi, width, ItrGaus, j_decay);
  }


}  // end namespace Chroma

