// $Id: unprec_slrc_fermact_w.cc,v 3.3 2007-03-05 20:03:16 bglaessle Exp $
/*! \file
 *  \brief Unpreconditioned SLRC fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/fermact_factory_w.h"

#include "actions/ferm/linop/unprec_slrc_linop_w.h"
#include "actions/ferm/fermacts/unprec_slrc_fermact_w.h"

//#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermstates/ferm_createstate_reader_w.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecSLRCFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,
		      multi1d<LatticeColorMatrix>,
		      multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
								     const std::string& path)
    {
      return new UnprecSLRCFermAct(CreateFermStateEnv::reader(xml_in, path), 
				     CloverFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion,
		  multi1d<LatticeColorMatrix>,
		  multi1d<LatticeColorMatrix> >* createFermAct(XMLReader& xml_in,
							       const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_SLRC";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
	success &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
	registered = true;
      }
      return success;
    }
  }


  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  UnprecLinearOperator<LatticeFermion,
		       multi1d<LatticeColorMatrix>,
		       multi1d<LatticeColorMatrix> >* 
  UnprecSLRCFermAct::linOp(Handle< FermState<T,P,Q> > state) const
  {
    return new UnprecSLRCLinOp(state,param);
  }

}

