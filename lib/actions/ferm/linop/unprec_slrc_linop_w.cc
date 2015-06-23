// $Id: unprec_slrc_linop_w.cc,v 3.0 2012-03-23 04:58:51 bglaessle Exp $
/*! \file
 *  \brief Unpreconditioned slrc linear operator
 */

#include "chromabase.h"
#include "actions/ferm/fermstates/simple_fermstate.h"
#include "actions/ferm/fermstates/periodic_fermstate.h"
#include "actions/ferm/fermstates/stout_fermstate_w.h"
#include "actions/ferm/linop/unprec_slrc_linop_w.h"

using namespace QDP::Hints;

namespace Chroma 
{ 
  //! Creation routine with Anisotropy
  /*!
   * \param fs 	    gauge field     	       (Read)
   * \param param_  parameters   	       (Read)
   */
  void UnprecSLRCLinOp::create(Handle< FermState<T,P,Q> > fs,
				 const CloverFermActParams& param_)
  {
    //   QDPIO::cout << __PRETTY_FUNCTION__ << ": enter" << endl;

    param = param_;

    // Need to make sure that fs is a stout ferm state
    // We want to have Clover with thin links
    thin_fs  = new PeriodicFermState<T,P,Q>( fs.cast< SLICFermState<T, P, Q> >()->getThinLinks() );

    A.create(thin_fs, param);
    D.create(fs, param.anisoParam);

    // QDPIO::cout << __PRETTY_FUNCTION__ << ": exit" << endl;
  }



  //! Apply unpreconditioned SLRC fermion linear operator
  /*!
   * The operator acts on the entire lattice
   *
   * \param chi 	  Pseudofermion field          (Write)
   * \param psi 	  Pseudofermion field          (Read)
   * \param isign   Flag ( PLUS | MINUS )   	       (Read)
   */
  void UnprecSLRCLinOp::operator()(LatticeFermion & chi, 
				     const LatticeFermion& psi, 
				     enum PlusMinus isign) const
  {
    LatticeFermion tmp; moveToFastMemoryHint(tmp);
    Real mhalf = -0.5;

    //  chi   =  A . psi - 0.5 * D' . psi  */
    A(chi, psi, isign); 
    D(tmp, psi, isign);
    chi += mhalf * tmp;

    getFermBC().modifyF(chi);
  }


  void 
  UnprecSLRCLinOp::deriv(multi1d<LatticeColorMatrix>& ds_u,
			   const LatticeFermion& chi, const LatticeFermion& psi, 
			   enum PlusMinus isign) const
  {
    // A. deriv will resize
    
    A.deriv(ds_u, chi, psi, isign);

    multi1d<LatticeColorMatrix> ds_tmp(Nd);

    ds_tmp = zero;
    D.deriv(ds_tmp, chi, psi, isign);
    for(int mu=0; mu < Nd; mu++) { 
      ds_u[mu] -= Real(0.5)*ds_tmp[mu];
    }
    
    getFermBC().zero(ds_u);
  }


  //! Return flops performed by the operator()
  unsigned long UnprecSLRCLinOp::nFlops() const
  {
    unsigned long site_flops = D.nFlops()+A.nFlops()+4*Nc*Ns;
    return site_flops*Layout::sitesOnNode();
  }

} // End Namespace Chroma
