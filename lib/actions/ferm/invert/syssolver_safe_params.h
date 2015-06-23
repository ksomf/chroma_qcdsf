// -*- C++ -*-
// $Id: syssolver_safe_params.h, 20012-03-06 15:20:54 bglaessle jnajjar Exp $
/*! \file
 *  \brief Solve a CG1 system
 */

#ifndef __syssolver_safe_params_h__
#define __syssolver_safe_params_h__

#include "chromabase.h"


namespace Chroma
{

  //! Params for CG inverter
  /*! \ingroup invert */
  struct SysSolverSafeParams
  {
    SysSolverSafeParams();
    SysSolverSafeParams(XMLReader& in, const std::string& path);

    Real          RsdCG;           /*!< CG residual */
    Real          RsdBiCGStab;           /*!< CG residual */
    int           MaxCG;           /*!< Maximum CG iterations */
    int           MaxBiCGStab;           /*!< Maximum CG iterations */
    int           MinCG;           /*!< Minimum CG iterations (useful for charm) */

    int 	  Strategy;

    Real          RsdCGRestart;    /*!< CG residual for a possibly double precision restart. Only valid for some solvers eg CG-DWF */

    int           MaxCGRestart;    /*!< Max no of CG iterations for a possibly double precision restart. Only valid for some solvers, eg CG-DWF */

    Real max_BiCG_relative_Res; /*!< In test for real convergence: Discard solves with higher relative residual than this */
  };


  // Reader/writers
  /*! \ingroup invert */
  void read(XMLReader& xml, const std::string& path, SysSolverSafeParams& param);

  /*! \ingroup invert */
  void write(XMLWriter& xml, const std::string& path, const SysSolverSafeParams& param);

} // End namespace

#endif
