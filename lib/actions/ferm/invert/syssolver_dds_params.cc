// $Id: syssolver_cg_params.cc,v 3.5 2008-03-25 10:43:44 mcneile Exp $
/*! \file
 *  \brief Params of CG inverter
 */

#include "actions/ferm/invert/syssolver_dds_params.h"

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const std::string& path, SysSolverDDSParams& param)
  {
    XMLReader paramtop(xml, path);

    param.BS.resize(Nd);

    read(paramtop, "RsdDDS", param.RsdDDS);
    read(paramtop, "MaxDDS", param.MaxDDS);
    read(paramtop, "Kappa", param.Kappa);
    read(paramtop, "Csw", param.Csw);

    read(paramtop, "Nkv", param.Nkv);
    read(paramtop, "Ncy", param.Ncy);
    read(paramtop, "Nmr", param.Nmr);
    read(paramtop, "DeflatedNV", param.DeflatedNV);
    read(paramtop, "BlkRel", param.BlkRel);
    read(paramtop, "BlockSize", param.BS);
    read(paramtop, "SLRC", param.slrc);
    /*
    Real          RsdDDS;           !<  residual
    Real          Kappa;           !<  kappa
    Real          Csw;           !< csw
    int           MaxDDS;           !< Maximum outer iterations
    int           Nkv; // max krylov space dimension
    int           Ncy; // SAP cycles
    int           Nmr; //# relaxation iterations on the blocks
    int           DeflatedNV; // 0 = fgcr, 1 = fgmres, >=2 = fgmres-dr with "DeflatedNV" deflated vectors
    int           BlkRel;   // 0 = MR, 1 = Gauss-Seidel, 2 = poly (not yet implemented)
    multi1d<int>  bs;
    */

   /* if( paramtop.count("RsdCGRestart") > 0 ) {
      read(paramtop, "RsdCGRestart", param.RsdCGRestart);
    }
    else {
      param.RsdCGRestart = param.RsdCG;
    } */

   /* if( paramtop.count("MaxCGRestart") > 0 ) {
      read(paramtop, "MaxCGRestart", param.MaxCGRestart);
    }
    else {
      param.MaxCGRestart = param.MaxCG;
    } */


    /* int aa = paramtop.count("MinCG") ;
    if( aa  > 0 ) {
      read(paramtop,"MinCG",param.MinCG);
    }
    else {
      param.MinCG = 0 ;
    }
    */

  }

  // Writer parameters
  void write(XMLWriter& xml, const std::string& path, const SysSolverDDSParams& param)
  {
    push(xml, path);

//  int version = 1;
//  write(xml, "version", version);
    write(xml, "invType", "DDS_INVERTER");
    write(xml, "RsdDDS", param.RsdDDS);
    write(xml, "MaxDDS", param.MaxDDS);
    write(xml, "Kappa", param.Kappa);
    write(xml, "Csw", param.Csw);

 //   write(xml, "MinCG", param.MinCG);
 //   write(xml, "RsdCGRestart", param.RsdCGRestart);
 //   write(xml, "MaxCGRestart", param.MaxCGRestart);
    pop(xml);
  }

  //! Default constructor
  SysSolverDDSParams::SysSolverDDSParams()
  {
    RsdDDS = zero;
    MaxDDS = 0;
    slrc = 0;  //default 0 = clover

  }

  //! Read parameters
  SysSolverDDSParams::SysSolverDDSParams(XMLReader& xml, const std::string& path)
  {
    read(xml, path, *this);
  }

}
