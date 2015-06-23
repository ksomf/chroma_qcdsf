#ifndef __SYSSOLVER_REL_BICGSTAB_PARAMS_SLRC_FEYNHELL_QCDSF_H__
#define __SYSSOLVER_REL_BICGSTAB_PARAMS_SLRC_FEYNHELL_QCDSF_H__

#include "chromabase.h"
#include "actions/ferm/fermacts/clover_fermact_params_w.h"
#include "actions/ferm/fermacts/slrc_feynhell_fermact_params_w.h"
#include "io/xml_group_reader.h"

namespace Chroma 
{
  struct SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF { 
    SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF(XMLReader& xml, const std::string& path);
    SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF() {};
    SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF( const SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF& p) {
      clovParams = p.clovParams;
      fhParams = p.fhParams;
      MaxIter = p.MaxIter;
      RsdTarget = p.RsdTarget;
      Delta = p.Delta;
    }
    CloverFermActParams clovParams;
    SLRCFeynHellFermActParams fhParams;
    int MaxIter;
    Real RsdTarget;
    Real Delta;
  };

  typedef SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF SysSolverReliableCGCloverParamsFeynHellQCDSF;
  void read(XMLReader& xml, const std::string& path, SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF& p);

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF& param);



}

#endif


