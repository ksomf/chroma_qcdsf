// -*- C++ -*-
// $Id: unprec_slrc_feynhell_dumb_linop_w.h,v 3.0 2012-03-23 04:58:51 bglaessle Exp $
/*! \file
 *  \brief Unpreconditioned SLRC fermion linear operator
 */

#ifndef __unprec_slrc_feynhell_dumb_linop_w_h__
#define __unprec_slrc_feynhell_dumb_linop_w_h__

#include "linearop.h"
#include "actions/ferm/linop/dslash_w.h"
#include "actions/ferm/linop/clover_term_w.h"
#include "actions/ferm/fermacts/slrc_feynhell_fermact_params_w.h"

namespace Chroma
{
  //! Unpreconditioned Slrc-Dirac operator
  /*!
   * \ingroup linop
   *
   * This routine is specific to Wilson fermions!
   */
  
  class UnprecSlrcFeynHellDumbFLinOpQCDSF : public UnprecLinearOperator<LatticeFermionF, 
			    multi1d<LatticeColorMatrixF>, multi1d<LatticeColorMatrixF> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermionF T;
    typedef LatticeColorMatrixF U;
    typedef multi1d<U>  P;
    typedef multi1d<U>  Q;

    //! Partial constructor
    UnprecSlrcFeynHellDumbFLinOpQCDSF() {}

    //! Full constructor
    UnprecSlrcFeynHellDumbFLinOpQCDSF(Handle< FermState<T,P,Q> > fs,
				      Handle< FermState<T,P,Q> > thin_fs,
			    const CloverFermActParams& param_,
			    const SLRCFeynHellFermActParams& fhparam_)
    {create(fs,thin_fs,param_,fhparam_);}
    
    //! Destructor is automatic
    ~UnprecSlrcFeynHellDumbFLinOpQCDSF() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
		Handle< FermState<T,P,Q> > thin_fs,
		const CloverFermActParams& param_,
		const SLRCFeynHellFermActParams& fhparam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermionF& chi, const LatticeFermionF& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned SLRC dM/dU
    void deriv(multi1d<LatticeColorMatrixF>& ds_u, 
	       const LatticeFermionF& chi, const LatticeFermionF& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    //    Handle< FermState<T,P,Q> > thin_fs;
    CloverFermActParams param;
    SLRCFeynHellFermActParams fhparam;
    WilsonDslashF        D;
    CloverTermF          A;
  };

  class UnprecSlrcFeynHellDumbDLinOpQCDSF : public UnprecLinearOperator<LatticeFermionD, 
			    multi1d<LatticeColorMatrixD>, multi1d<LatticeColorMatrixD> >
  {
  public:
    // Typedefs to save typing
    typedef LatticeFermionD T;
    typedef LatticeColorMatrixD U;
    typedef multi1d<U>  P;
    typedef multi1d<U>  Q;

    //! Partial constructor
    UnprecSlrcFeynHellDumbDLinOpQCDSF() {}

    //! Full constructor
    UnprecSlrcFeynHellDumbDLinOpQCDSF(Handle< FermState<T,P,Q> > fs,
				      Handle< FermState<T,P,Q> > thin_fs,
			    const CloverFermActParams& param_,
			    const SLRCFeynHellFermActParams& fhparam_)
    {create(fs,thin_fs,param_,fhparam_);}
    
    //! Destructor is automatic
    ~UnprecSlrcFeynHellDumbDLinOpQCDSF() {}

    //! Return the fermion BC object for this linear operator
    const FermBC<T,P,Q>& getFermBC() const {return D.getFermBC();}

    //! Creation routine
    void create(Handle< FermState<T,P,Q> > fs,
	      Handle< FermState<T,P,Q> > thin_fs,
		const CloverFermActParams& param_,
		const SLRCFeynHellFermActParams& fhparam_);

    //! Apply the operator onto a source vector
    void operator() (LatticeFermionD& chi, const LatticeFermionD& psi, enum PlusMinus isign) const;

    //! Derivative of unpreconditioned SLRC dM/dU
    void deriv(multi1d<LatticeColorMatrixD>& ds_u, 
	       const LatticeFermionD& chi, const LatticeFermionD& psi, 
	       enum PlusMinus isign) const;

    //! Return flops performed by the operator()
    unsigned long nFlops() const;

  private:
    //    Handle< FermState<T,P,Q> > thin_fs;
    CloverFermActParams param;
    SLRCFeynHellFermActParams fhparam;
    WilsonDslashD        D;
    CloverTermD          A;
  };


} // End Namespace Chroma


#endif
