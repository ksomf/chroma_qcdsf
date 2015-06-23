// $Id: syssolver_safe_params.cc, 20012-03-06 15:20:54 bglaessle jnajjar Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_safe_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, SysSolverSafeParams& param)
  {
    XMLReader paramtop(xml, path);

    read(paramtop, "RsdBiCGStab", param.RsdBiCGStab);
    read(paramtop, "MaxBiCGStab", param.MaxBiCGStab);
    read(paramtop, "RsdCG", param.RsdCG);
    read(paramtop, "MaxCG", param.MaxCG);
    read(paramtop, "Strategy", param.Strategy);

    if( paramtop.count("RsdCGRestart") > 0 ) { 
      read(paramtop, "RsdCGRestart", param.RsdCGRestart);
    }
    else {
      param.RsdCGRestart = param.RsdCG;
    }

    if( paramtop.count("MaxCGRestart") > 0 ) { 
      read(paramtop, "MaxCGRestart", param.MaxCGRestart);
    }
    else {
      param.MaxCGRestart = param.MaxCG;
    }


    int aa = paramtop.count("MinCG") ;
    if( aa  > 0 ) { 
      read(paramtop,"MinCG",param.MinCG);
    }
    else {
      param.MinCG = 0 ; 
    }
    
    if( paramtop.count("max_BiCG_relative_Res") > 0 ) { 
    read(paramtop, "max_BiCG_relative_Res",param.max_BiCG_relative_Res);
    }
    else param.max_BiCG_relative_Res=1e-4;
  }

  // Writer parameters
  void write(XMLWriter& xml, const string& path, const SysSolverSafeParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "SAFE_INVERTER");
    write(xml, "RsdBiCGStab", param.RsdBiCGStab);
    write(xml, "MaxBiCGStab", param.MaxBiCGStab);
    write(xml, "RsdCG", param.RsdCG);
    write(xml, "MaxCG", param.MaxCG);
    write(xml, "MinCG", param.MinCG);
    write(xml, "RsdCGRestart", param.RsdCGRestart);
    write(xml, "MaxCGRestart", param.MaxCGRestart);
    write(xml, "Strategy", param.Strategy);
    pop(xml);
  }

  //! Default constructor
  SysSolverSafeParams::SysSolverSafeParams()
  {
    RsdCG = zero;
    RsdBiCGStab = zero;
    MaxCG = 0;
    MinCG = 0;
    RsdCGRestart = RsdCG;
    MaxCGRestart = MaxCG;
  }

  //! Read parameters
  SysSolverSafeParams::SysSolverSafeParams(XMLReader& xml, const string& path)
  {
    read(xml, path, *this);
  }

}
