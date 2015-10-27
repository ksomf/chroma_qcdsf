#include "actions/ferm/invert/syssolver_rel_bicgstab_slrc_feynhell_params_qcdsf.h"
#include "chromabase.h"
#include "io/xml_group_reader.h"



using namespace QDP;

namespace Chroma {

  SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF::SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF(XMLReader& xml,
						       const std::string& path)
  {
    QDPIO::cout << "Creating SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF" << std::endl;
    QDPIO::cout << "XML path is " << path << std::endl;

    XMLReader paramtop(xml, path);
    read(paramtop, "MaxIter", MaxIter);
    read(paramtop, "RsdTarget", RsdTarget);
    read(paramtop, "CloverParams", clovParams);
    read(paramtop, path, fhParams);
    read(paramtop, "Delta", Delta);
  }

  void read(XMLReader& xml, const std::string& path,
	    SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF& p)
  {
    QDPIO::cout << "Reading SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF" << std::endl;

    SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF tmp(xml, path);
    p = tmp;
  }

  void write(XMLWriter& xml, const std::string& path,
	     const SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF& p) {

    QDPIO::cout << "Writing SysSolverReliableBiCGStabSlrcFeynHellParamsQCDSF" << std::endl;
    push(xml, path);
    write(xml, "MaxIter", p.MaxIter);
    write(xml, "RsdTarget", p.RsdTarget);
    write(xml, "CloverParams", p.clovParams);
    write(xml, "FeynHellParam", p.fhParams);
    write(xml, "Delta", p.Delta);
    pop(xml);

  }



}
