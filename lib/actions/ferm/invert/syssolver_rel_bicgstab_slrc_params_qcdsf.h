#ifndef __SYSSOLVER_REL_BICGSTAB_PARAMS_SLRC_QCDSF_H__
#define __SYSSOLVER_REL_BICGSTAB_PARAMS_SLRC_QCDSF_H__

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{
  struct SysSolverReliableBiCGStabSlrcParamsQCDSF { 
    SysSolverReliableBiCGStabSlrcParamsQCDSF(XMLReader& xml, const std::string& path);
    SysSolverReliableBiCGStabSlrcParamsQCDSF() {};
    SysSolverReliableBiCGStabSlrcParamsQCDSF( const SysSolverReliableBiCGStabSlrcParamsQCDSF& p) {
      clovParams = p.clovParams;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
    }
    CloverFermActParams clovParams;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
  };

  typedef SysSolverReliableBiCGStabSlrcParamsQCDSF SysSolverReliableCGCloverParamsQCDSF;
  void read(XMLReader& xml, const std::string& path, SysSolverReliableBiCGStabSlrcParamsQCDSF& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverReliableBiCGStabSlrcParamsQCDSF& param);



}

#endif


