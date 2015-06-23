// $Id: eoprec_slrc_feynhell_linop_w.cc,v 3.8 2007-11-29 14:19:41 bjoo Exp $
/*! \file
 *  \brief Even-odd preconditioned Clover linear operator (fat-relevant, thin-irrelevant terms)
 *
 * Here, the relevant terms are smeared and the irrelevant terms are not smeared.
 * Code provided by Thomas Kaltenbrunner.
 *
 */

#include "actions/ferm/linop/eoprec_slrc_feynhell_linop_w.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
using namespace QDP::Hints;

namespace Chroma
{

  //! Creation routine with Anisotropy
  /*!
   * \param u_ 	    gauge field     	       (Read)
   * \param param_  fermion kappa   	       (Read)
   */
  void EvenOddPrecSLRCFeynHellLinOp::create(Handle< FermState<T,P,Q> > fs,
					    const CloverFermActParams& param_,
					    const SLRCFeynHellFermActParams& fhparam_)
  {
    START_CODE();

    // Alex check
    QDPIO::cout << "I'm the linear operator Alex made!!" << std::endl;

    param = param_;
    slrc_fs = fs; // added by I.Kanamori
    fhparam = fhparam_;
    QDPIO::cout << "(debug info) a local version of create is called!\n";

    // Need to make sure that fs is a stout ferm state
    // We want to have Clover with thin links
    thin_fs  = new PeriodicFermState<T,P,Q>(
					    	fs.cast< SLICFermState<T, P, Q> >()->getThinLinks());

    clov.create(thin_fs, param);

    invclov.create(thin_fs,param,clov);  // make a copy
    invclov.choles(0);  // invert the cb=0 part

    //...and WilsonDslash with stout links
    D.create(fs, param.anisoParam);

    END_CODE();
  }

  //! Apply the the odd-odd block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::oddOddLinOp(LatticeFermion& chi, const LatticeFermion& psi,
				      enum PlusMinus isign) const
  {
    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
	Complex phalf = 0.5;// * fhparam.lambda;

    clov.apply(chi, psi, isign, 1);

    //tmp1[rb[1]] = Gamma(fhparam.op) * psi;
	tmp1[rb[1]] = psi;

    if( isign == PLUS ) {
      chi[rb[1]] += phalf * tmp1;
    }
    else {
      chi[rb[1]] -= phalf * tmp1;
    }
  }


  //! Apply the the even-even block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::evenEvenLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					enum PlusMinus isign) const
  {
    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
	Complex phalf = 0.5;// * fhparam.lambda;

    // Nuke for testing
    clov.apply(chi, psi, isign, 0);

    //tmp1[rb[0]] = Gamma(fhparam.op) * psi;
	tmp1[rb[0]] = psi;

    if( isign == PLUS ) {
      chi[rb[0]] += phalf * tmp1;
    }
    else {
      chi[rb[0]] -= phalf * tmp1;
    }
  }

  //! Apply the inverse of the even-even block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::evenEvenInvLinOp(LatticeFermion& chi, const LatticeFermion& psi,
					   enum PlusMinus isign) const
  {
    invclov.apply(chi, psi, isign, 0);
  }


  //! Apply even-odd linop component
  /*!
   * The operator acts on the entire even sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  EvenOddPrecSLRCFeynHellLinOp::evenOddLinOp(LatticeFermion& chi,
				       const LatticeFermion& psi,
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 0);
    chi[rb[0]] *= mhalf;

    END_CODE();
  }

  //! Apply odd-even linop component
  /*!
   * The operator acts on the entire odd sublattice
   *
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void
  EvenOddPrecSLRCFeynHellLinOp::oddEvenLinOp(LatticeFermion& chi,
				       const LatticeFermion& psi,
				       enum PlusMinus isign) const
  {
    START_CODE();

    Real mhalf = -0.5;

    D.apply(chi, psi, isign, 1);
    chi[rb[1]] *= mhalf;

    END_CODE();
  }


  //! Apply even-odd preconditioned Clover fermion linear operator
  /*!
   * \param chi 	  Pseudofermion field     	       (Write)
   * \param psi 	  Pseudofermion field     	       (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void EvenOddPrecSLRCFeynHellLinOp::operator()(LatticeFermion & chi,
					  const LatticeFermion& psi,
					  enum PlusMinus isign) const
  {
    START_CODE();

    LatticeFermion tmp1; moveToFastMemoryHint(tmp1);
    LatticeFermion tmp2; moveToFastMemoryHint(tmp2);
    Real mquarter = -0.25;
	Complex phalf = 0.5; //* fhparam.lambda;

    //  tmp1_o  =  D_oe   A^(-1)_ee  D_eo  psi_o
    D.apply(tmp1, psi, isign, 0);
    invclov.apply(tmp2, tmp1, isign, 0);
    D.apply(tmp1, tmp2, isign, 1);

    //  chi_o  =  A_oo  psi_o  -  tmp1_o
    clov.apply(chi, psi, isign, 1);

    chi[rb[1]] += mquarter*tmp1;

    //    tmp1[rb[1]] = Gamma(fhparam.op) * psi;
    //tmp1 = Gamma(fhparam.op) * psi;
		tmp1[rb[1]] = psi;
	tmp1 = psi;

    //    if( isign == PLUS ) {
    chi[rb[1]] += phalf * tmp1;
    //    }
    //    else {
    //      chi[rb[1]] -= phalf * tmp1;
    //    }

    END_CODE();
  }


  //! Apply the even-even block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::derivEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
					     const LatticeFermion& chi, const LatticeFermion& psi,
					     enum PlusMinus isign) const
  {
    START_CODE();

    clov.deriv(ds_u, chi, psi, isign, 0);
    getFermBC().zero(ds_u);

    END_CODE();
  }

  //! Apply the even-even block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::derivLogDetEvenEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
						 enum PlusMinus isign) const
  {
    START_CODE();

    invclov.derivTrLn(ds_u, isign, 0);
    getFermBC().zero(ds_u);

    END_CODE();
  }

  //! Apply the the even-odd block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::derivEvenOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, const LatticeFermion& psi,
					    enum PlusMinus isign) const
  {
    START_CODE();

    //temp variable for the fat force
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    // Dslash will resize this.
    ds_u.resize(Nd);


    D.deriv(ds_tmp, chi, psi, isign, 0);
    for(int mu=0; mu < Nd; mu++) {
      ds_tmp[mu]  *= Real(-0.5);
    }

    //WilsonDslash is thick, so need fatForceToThin here! (includes change of BCs)
    SLICFermState<T, P, Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slrc_fs);
    sfs.fatForceToThin(ds_tmp,ds_u);
    getFermBC().zero(ds_u);

    END_CODE();
  }

  //! Apply the the odd-even block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::derivOddEvenLinOp(multi1d<LatticeColorMatrix>& ds_u,
					    const LatticeFermion& chi, const LatticeFermion& psi,
					    enum PlusMinus isign) const
  {
    START_CODE();

    //temp variable for the fat force
    multi1d<LatticeColorMatrix> ds_tmp(Nd);
    ds_u.resize(Nd);

    D.deriv(ds_tmp, chi, psi, isign, 1);
    for(int mu=0; mu < Nd; mu++) {
      ds_tmp[mu]  *= Real(-0.5);
    }

    SLICFermState<T, P, Q>& sfs = dynamic_cast<SLICFermState<T,P,Q>& >(*slrc_fs);
    sfs.fatForceToThin(ds_tmp,ds_u);
    getFermBC().zero(ds_u);

    END_CODE();
  }

  // Inherit this
  //! Apply the the odd-odd block onto a source vector
  void
  EvenOddPrecSLRCFeynHellLinOp::derivOddOddLinOp(multi1d<LatticeColorMatrix>& ds_u,
					   const LatticeFermion& chi, const LatticeFermion& psi,
					   enum PlusMinus isign) const
  {
    START_CODE();

    ds_u.resize(Nd);
    clov.deriv(ds_u, chi, psi, isign, 1);
    getFermBC().zero(ds_u);

    END_CODE();
  }


  //! Return flops performed by the operator()
  unsigned long EvenOddPrecSLRCFeynHellLinOp::nFlops() const
  {
    unsigned long cbsite_flops = 2*D.nFlops()+2*clov.nFlops()+4*Nc*Ns;
    return cbsite_flops*(Layout::sitesOnNode()/2);
  }

  //! Get the log det of the even even part
  // BUt for now, return zero for testing.
  Double EvenOddPrecSLRCFeynHellLinOp::logDetEvenEvenLinOp(void) const  {
    return invclov.cholesDet(0);
  }
} // End Namespace Chroma
