// $Id: syssolver_linop_rel_bicgstab_clover.cc,v 3.1 2009-05-20 15:25:52 bjoo Exp $
/*! \file
 *  \brief Solve a MdagM*psi=chi linear system by CG2
 */

#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"
#include "actions/ferm/invert/syssolver_rel_bicgstab_slrc_params_qcdsf.h"
#include "actions/ferm/invert/syssolver_linop_rel_bicgstab_slrc_qcdsf.h"

namespace Chroma
{
  namespace LinOpSysSolverReliableBiCGStabSlrcEnvQCDSF
  {

    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("RELIABLE_BICGSTAB_MP_SLRC_INVERTER-QCDSF");

      //! Local registration flag
      bool registered = false;
    }



    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,	
						  const std::string& path,
						  Handle< FermState< LatticeFermion, multi1d<LatticeColorMatrix>, multi1d<LatticeColorMatrix> > > state, 
						  
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverReliableBiCGStabSlrcQCDSF(A, state,SysSolverReliableBiCGStabSlrcParamsQCDSF(xml_in, path));
    }

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheLinOpFermSystemSolverFactory::Instance().registerObject(name, createFerm);
	registered = true;
      }
      return success;
    }
  }
}
