// -*- C++ -*-
// $Id: wuppertal_smear.h, v 1.0 2011-04-08 ehmann, rwschiel Exp $
/*! \file
 *  \brief Wuppertal smearing of color vector
 */

#ifndef __wuppertal_smear_h__
#define __wuppertal_smear_h__

namespace Chroma 
{

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
		       const Real& width, int ItrGaus, int j_decay);


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
		       const Real& width, int ItrGaus, int j_decay);


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
		       const Real& width, int ItrGaus, int j_decay);


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
		       const Real& width, int ItrGaus, int j_decay);

}  // end namespace Chroma

#endif
