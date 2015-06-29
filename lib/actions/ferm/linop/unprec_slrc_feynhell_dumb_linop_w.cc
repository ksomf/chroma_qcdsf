// $Id: unprec_slrc_feynhell_dumb_linop_w.cc,v 3.0 2012-03-23 04:58:51 bglaessle Exp $
/*! \file
 *  \brief Unpreconditioned slrc linear operator
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/linop/unprec_slrc_feynhell_dumb_linop_w.h"

#include "util/ft/sftmom.h"

using namespace QDP::Hints;

namespace Chroma
{
  //! Creation routine with Anisotropy
  /*!
   * \param fs 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecSlrcFeynHellDumbFLinOpQCDSF::create(Handle< FermState<T,P,Q> > fs,
						 Handle< FermState<T,P,Q> > thin_fs,
				       const CloverFermActParams& param_,
				       const SLRCFeynHellFermActParams& fhparam_)
  {
    QDPIO::cout << "Creating UnprecSlrcFeynHellDumbFLinOpQCDSF" << std::endl;

    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;
    fhparam = fhparam_;

    // Need to make sure that fs is a stout ferm state
    // We want to have Clover with thin links
//    thin_fs  = new PeriodicFermState<T,P,Q>( fs.cast< SLICFermState<T, P, Q> >()->getThinLinks() );

    A.create(thin_fs, param);
    D.create(fs, param.anisoParam);

    QDPIO::cout << "Done creating UnprecSlrcFeynHellDumbFLinOpQCDSF" << std::endl;

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
  }



  //! Apply unpreconditioned SLRC fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecSlrcFeynHellDumbFLinOpQCDSF::operator()(LatticeFermionF & chi,
				     const LatticeFermionF& psi,
				     enum PlusMinus isign) const
  {
    LatticeFermionF tmp; moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;
    Real phalf = 0.5;

    A(chi, psi, isign);
    D(tmp, psi, isign);
    chi += mhalf * tmp;

    for (int i = 0; i < fhparam.numParam; i++)
	{
		chi += phalf
			* fhparam.FHparam[i].lambda
			* fhparam.FHparam[i].noise
			* fhparam.FHparam[i].phases
			* (Gamma(fhparam.FHparam[i].op) * psi);
	}

    getFermBC().modifyF(chi);
  }


  void
  UnprecSlrcFeynHellDumbFLinOpQCDSF::deriv(multi1d<LatticeColorMatrixF>& ds_u,
			   const LatticeFermionF& chi, const LatticeFermionF& psi,
			   enum PlusMinus isign) const
  {
    // A. deriv will resize

    A.deriv(ds_u, chi, psi, isign);

    multi1d<LatticeColorMatrixF> ds_tmp(Nd);

    ds_tmp = zero;
    D.deriv(ds_tmp, chi, psi, isign);
    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] -= Real(0.5)*ds_tmp[mu];
    }

    getFermBC().zero(ds_u);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecSlrcFeynHellDumbFLinOpQCDSF::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+A.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }


  void UnprecSlrcFeynHellDumbDLinOpQCDSF::create(Handle< FermState<T,P,Q> > fs,
						 Handle< FermState<T,P,Q> > thin_fs,
				       const CloverFermActParams& param_,
				       const SLRCFeynHellFermActParams& fhparam_)
  {
    QDPIO::cout << "Creating UnprecSlrcFeynHellDumbDLinOpQCDSF" << std::endl;

    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << std::endl;

    param = param_;
    fhparam = fhparam_;

    // Need to make sure that fs is a stout ferm state
    // We want to have Clover with thin links
//    thin_fs  = new PeriodicFermState<T,P,Q>( fs.cast< SLICFermState<T, P, Q> >()->getThinLinks() );

    A.create(thin_fs, param);
    D.create(fs, param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << std::endl;
  }



  //! Apply unpreconditioned SLRC fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecSlrcFeynHellDumbDLinOpQCDSF::operator()(LatticeFermionD & chi,
				     const LatticeFermionD& psi,
				     enum PlusMinus isign) const
  {
    LatticeFermionD tmp; moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;
    Real phalf = 0.5;

    A(chi, psi, isign);
    D(tmp, psi, isign);
    chi += mhalf * tmp;

    for (int i = 0; i < fhparam.numParam; i++)
	{
		chi += phalf
			* fhparam.FHparam[i].lambda
			* fhparam.FHparam[i].noise
			* fhparam.FHparam[i].phases
			* (Gamma(fhparam.FHparam[i].op) * psi);
	}

    getFermBC().modifyF(chi);
  }


  void
  UnprecSlrcFeynHellDumbDLinOpQCDSF::deriv(multi1d<LatticeColorMatrixD>& ds_u,
			   const LatticeFermionD& chi, const LatticeFermionD& psi,
			   enum PlusMinus isign) const
  {
    // A. deriv will resize

    A.deriv(ds_u, chi, psi, isign);

    multi1d<LatticeColorMatrixD> ds_tmp(Nd);

    ds_tmp = zero;
    D.deriv(ds_tmp, chi, psi, isign);
    for(int mu=0; mu < Nd; mu++) {
      ds_u[mu] -= Real(0.5)*ds_tmp[mu];
    }

    getFermBC().zero(ds_u);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecSlrcFeynHellDumbDLinOpQCDSF::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+A.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
