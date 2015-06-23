// -*- C++ -*-
// $Id: wuppertal_source_smear.h, v 1.0 2011-08-31 15:58:33 bglaessle Exp $
/*! \file
 *  \brief Wuppertal smearing for the SOURCE:
 		- only a subset (=timeslice) of the lattice is smeared
 		- routine for spin-diagonal SOURCE: only one spin component of 16 is smeared an reinserted
 */

#ifndef __wuppertal_source_smear_h__
#define __wuppertal_source_smear_h__

namespace Chroma 
{
//  extern Set TSliceSetQCDSF;
//  extern Set TSliceSetQCDSF;
//  void initTimeSliceSet();

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
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice );


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
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice );


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
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice );


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
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice );

  void wuppertalSourceSmear_spinDiagonal (const multi1d<LatticeColorMatrix>& u, 
		       LatticePropagator& chi, 
		       const Real& width, int ItrGaus, int j_decay, const Subset& tslice );
		       

}  // end namespace Chroma

#endif
