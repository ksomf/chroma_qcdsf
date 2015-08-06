// -*- C++ -*-

//  syssolver_dds_params.h,v  Andrea Nobile, Luca Castagnini
/*! \file
 *
 */

#ifndef __syssolver_dds_params_h__
#define __syssolver_dds_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for CG inverter
  /*! \ingroup invert */
  struct SysSolverDDSParams
  {
    SysSolverDDSParams();
    SysSolverDDSParams(XMLReader& in, const std::string& path);

    int           slrc;              /* 0 = normal clover,  1 = SLRC */
    Real          RsdDDS;           /*!<  residual */
    Real          Kappa;           /*!<  kappa */
    Real          Csw;           /*!< csw */
    int           MaxDDS;           /*!< Maximum outer iterations */
    int           Nkv; // max krylov space dimension
    int           Ncy; // SAP cycles
    int           Nmr; //# relaxation iterations on the blocks
    int           DeflatedNV; // 0 = fgcr, 1 = fgmres, >=2 = fgmres-dr with "DeflatedNV" deflated vectors
    int           BlkRel;   // 0 = MR, 1 = Gauss-Seidel, 2 = poly (not yet implemented)
    multi1d<int>  BS;
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const std::string& path, SysSolverDDSParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const SysSolverDDSParams& param);

} // End namespace

#endif
