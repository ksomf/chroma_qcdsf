//  syssolver_linop_dds.cc,v Andrea Nobile, Luca Castagnini 
/*! \file
 *  \brief Solve a M*psi=chi linear system by DDS
 */
#include "state.h"
#include "actions/ferm/invert/syssolver_linop_factory.h"
#include "actions/ferm/invert/syssolver_linop_aggregate.h"

#include "actions/ferm/invert/syssolver_linop_dds.h"

namespace Chroma
{

  //! CG1 system solver namespace
  namespace LinOpSysSolverDDSEnv
  {
    //! Anonymous namespace
    namespace
    {
      //! Name to be used
      const std::string name("DDS_INVERTER");

      //! Local registration flag
      bool registered = false;
    }


    //! Callback function
    LinOpSystemSolver<LatticeFermion>* createFerm(XMLReader& xml_in,
						  const std::string& path,
						  Handle< FermState<
						  LatticeFermion, 
						  multi1d<LatticeColorMatrix>,
						  multi1d<LatticeColorMatrix> 
						  > 
						  > state, 
      
						  Handle< LinearOperator<LatticeFermion> > A)
    {
      return new LinOpSysSolverDDS<LatticeFermion>(A, state, SysSolverDDSParams(xml_in, path));
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
