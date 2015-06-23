#include "actions/ferm/invert/syssolver_rel_bicgstab_slrc_params_qcdsf.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {
  
  SysSolverReliableBiCGStabSlrcParamsQCDSF::SysSolverReliableBiCGStabSlrcParamsQCDSF(XMLReader& xml, 
						       const std::string& path)
  {
    XMLReader paramtop(xml, path);
    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", clovParams);
    read(paramtop, "Delta", Delta);
  }

  void read(XMLReader& xml, const std::string& path, 
	    SysSolverReliableBiCGStabSlrcParamsQCDSF& p)
  {
    SysSolverReliableBiCGStabSlrcParamsQCDSF tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path, 
	     const SysSolverReliableBiCGStabSlrcParamsQCDSF& p) {
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.clovParams);
    write(xml, "Delta", p.Delta);
    pop(xml);

  }



}
