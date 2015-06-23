// $Id: wuppertal_source_smear.cc, v 1.0 2011-08-31 15:58:33 bglaessle Exp $
/*! \file
 *  \brief Wuppertal smearing for the SOURCE:
 		- only a subset (=timeslice) of the lattice is smeared
 		- routine for spin-diagonal SOURCE: only one spin component of 16 is smeared an reinserted
 */

#include "chromabase.h"
#include "meas/smear/wuppertal_source_smear.h"
#include "util/ferm/transf.h"
//#include "actions/boson/operator/klein_gord.h"

namespace Chroma 
{
/*  */

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
  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       T& chi, 
		       const Real& kappa, int ItrJacobi, int j_decay, const Subset& tslice )
  { 
    T psi;
    Real norm;

     if (j_decay < Nd)
      norm = Real (1) + Real (2 * Nd - 2) * kappa;
     else
      norm = Real (1) + Real (2 * Nd) * kappa;
    
    for (int n = 0; n < ItrJacobi; n++)
    {
      psi[tslice] = chi;

     for (int mu = 0; mu < Nd; mu++)
      if (mu != j_decay)
      {
	chi[tslice] += kappa * (u [mu] * shift (psi, FORWARD, mu) + shift (adj (u [mu]) * psi, BACKWARD, mu));
      }
      chi[tslice] = chi / norm;
    } 
  }

/*  template<typename T>
  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       T& chi, 
		       const Real& kappa, int ItrJacobi, int j_decay, const Subset& tslice )
  { 
    T psi, tmp;
    Real norm;

     if (j_decay < Nd)
      norm = Real (1) + Real (2 * Nd - 2) * kappa;
     else
      norm = Real (1) + Real (2 * Nd) * kappa;
    
    for (int n = 0; n < ItrJacobi; n++)
    {
      psi[tslice] = chi;

     for (int mu = 0; mu < Nd; mu++)
      if (mu != j_decay)
      {
//		tmp[tslice] = adj (u [mu]) * psi;
//		chi[tslice] += kappa * (u [mu] * shift (psi, FORWARD, mu) + shift (tmp, BACKWARD, mu));

		tmp[tslice] = u [mu] * shift (psi, FORWARD, mu);
		chi[tslice] += kappa * (tmp + shift (adj (u [mu]) * psi, BACKWARD, mu));
      }
      chi[tslice] = chi / norm;
    } 
  } */

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

  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeColorVector& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice )
  {
    wuppertalSourceSmear<LatticeColorVector> (u, chi, width, ItrGaus, j_decay, tslice);
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

  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeFermion& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice )
  {
    wuppertalSourceSmear<LatticeFermion> (u, chi, width, ItrGaus, j_decay, tslice);
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

  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticeStaggeredPropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice )
  {
    wuppertalSourceSmear<LatticeStaggeredPropagator> (u, chi, width, ItrGaus, j_decay, tslice);
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

  void wuppertalSourceSmear (const multi1d<LatticeColorMatrix>& u, 
		       LatticePropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice )
  {
	wuppertalSourceSmear<LatticePropagator> (u, chi, width, ItrGaus, j_decay, tslice);
  }

  void wuppertalSourceSmear_spinDiagonal (const multi1d<LatticeColorMatrix>& u,
		       LatticePropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice )
  {
	int spin_in = 0; // can be 0 <= any < Ns !!
  	LatticeFermion ferm;
	LatticeColorVector cv;

	for( int color=0; color<Nc; ++color ) {
		PropToFerm( chi, ferm, color, spin_in );
		cv = peekSpin( ferm, spin_in );
		wuppertalSourceSmear<LatticeColorVector> (u, cv, width, ItrGaus, j_decay, tslice);
		
		for( int spin_out=0; spin_out<Ns; ++spin_out ){
			ferm = zero;
			pokeSpin( ferm, cv, spin_out );
			FermToProp( ferm, chi, color, spin_out );
		}
	}
  }


}  // end namespace Chroma

