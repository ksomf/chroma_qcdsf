// -*- C++ -*-
// $Id: slrc_feynhell_fermact_params_w.h,v 3.3 2008-11-10 17:59:07 bjoo Exp $
/*! \file
 *  \brief Parameters for Clover fermion action
 */

#ifndef __slrc_feynhell_fermact_params_w_h__
#define __slrc_feynhell_fermact_params_w_h__

#include "io/aniso_io.h"
#include "util/ft/sftmom.h"

namespace Chroma
{
  //! Params for clover ferm acts
  /*! \ingroup fermacts */
  struct SLRCFeynHellFermActParams
  {
    SLRCFeynHellFermActParams();
    SLRCFeynHellFermActParams(XMLReader& in, const std::string& path);

    int numParam;

    struct FHParam
    {
      Complex lambda;
      int op;
      multi1d<int> mom;
      multi1d<int> source;
      Complex noise;
      LatticeComplex phases;
    };
    multi1d<FHParam> FHparam;
  };


  // Reader/writers
  /*! \ingroup fermacts */
  void read(XMLReader& xml, const std::string& path, SLRCFeynHellFermActParams& param);
  void read(XMLReader& xml, const std::string& path, SLRCFeynHellFermActParams::FHParam& param);

  /*! \ingroup fermacts */
  void write(XMLWriter& xml, const std::string& path, const SLRCFeynHellFermActParams& param);
  void write(XMLWriter& xml, const std::string& path, const SLRCFeynHellFermActParams::FHParam& param);
}

#endif
