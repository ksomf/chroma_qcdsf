// -*- C++ -*-
// $Id: wuppertal_quark_smearing.h, v 1.0 2011-04-08 ehmann, rwschiel $
/*! \file
 *  \brief Wuppertal smearing of color vector and propagator
 */

#ifndef __wuppertal_quark_smearing_h__
#define __wuppertal_quark_smearing_h__

#include "meas/smear/quark_smearing.h"

namespace Chroma
{
  //! Name and registration
  /*! @ingroup smear */
  namespace WuppertalQuarkSmearingEnv
  {
    extern const std::string name;
    bool registerAll ();


    //! Params for Wuppertal quark smearing
    /*! @ingroup smear */
    struct Params
    {
      Params () {}
      Params (XMLReader& in, const std::string& path);
      void writeXML (XMLWriter& in, const std::string& path) const;

      Real wvf_param;                   /*!< Smearing width */
      int  wvfIntPar;                   /*!< Number of smearing hits */
      int  no_smear_dir;		/*!< No smearing in this direction */
    };


    //! Wuppertal quark smearing
    /*! @ingroup smear
     *
     * Wuppertal quark smearing object
     */
    template<typename T>
    class QuarkSmear : public QuarkSmearing<T>
    {
    public:
      //! Full constructor
      QuarkSmear (const Params& p) : params (p) {}

      //! Smear the quark
      void operator() (T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      QuarkSmear () {}

    private:
      Params  params;   /*!< smearing params */
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read (XMLReader& xml, const std::string& path, WuppertalQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write (XMLWriter& xml, const std::string& path, const WuppertalQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
