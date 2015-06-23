// -*- C++ -*-
// $Id: unprec_slrc_linop_w.h,v 3.0 2012-03-23 04:58:51 bglaessle Exp $
/*! \file
 *  \brief Unpreconditioned SLRC fermion linear operator
 */

#ifndef __unprec_slrc_linop_w_h__
#define __unprec_slrc_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"


namespace Chroma 
{ 
  //! Unpreconditioned SLRC-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  
  class UnprecSLRCLinOp : public UnprecLinearOperator<LatticeFermion, 
			    multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermion               T;
    typedef multi1d<LatticeColorMatrix>  P;
    typedef multi1d<LatticeColorMatrix>  Q;

    //! Partial constructor
    UnprecSLRCLinOp() {}

    //! Full constructor
    UnprecSLRCLinOp(Handle< FermState<T,P,Q> > fs,
		      const CloverFermActParams& param_)
      {create(fs,param_);}
    
    //! Destructor is automatic
    ~UnprecSLRCLinOp() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		const CloverFermActParams& param_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermion& chi, const LatticeFermion& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned SLRC dM/dU
    void deriv(multi1d<LatticeColorMatrix>& ds_u, 
	       const LatticeFermion& chi, const LatticeFermion& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    Handle< FermState<T,P,Q> > thin_fs;
    CloverFermActParams param;
    WilsonDslash        D;
    CloverTerm          A;
  };

} // End Namespace Chroma


#endif
