// -*- C++ -*-
/* \file
 * \brief Stout params
 */

#ifndef _stout_fermstate_params_h_
#define _stout_fermstate_params_h_

#include "chromabase.h"
#include "actions/ferm/fermacts/slrc_feynhell_fermact_params_w.h"

namespace Chroma 
{
  //! Params for stout-links
  /*! \ingroup fermstates */
  struct StoutFermStateParams
  {
    //! Default constructor
    StoutFermStateParams();
    StoutFermStateParams(XMLReader& in, const std::string& path);

    multi2d<Real>  rho;
    multi1d<bool>  smear_in_this_dirP;

    int            n_smear;
    Real           sm_fact;

    bool                      doing_fh;
    SLRCFeynHellFermActParams fh_params;
  };

  void read(XMLReader& xml, const std::string& path, StoutFermStateParams& p);
  void write(XMLWriter& xml, const std::string& path, const StoutFermStateParams& p);
}

#endif
