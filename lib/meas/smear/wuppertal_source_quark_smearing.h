// -*- C++ -*-
// $Id: wuppertal_source_quark_smearing.h, v 1.0 2011-08-31 16:19:59 bglaessle $
/*! \file
 *  \brief \brief Wuppertal SOURCE smearing "operator", optimizes:
 		- if source(propagator) is spin-diagonal, only 1 of 16 spin components is smeared
 		- if source is on a single timeslice, only that is smeared
 */

#ifndef __wuppertal_source_quark_smearing_h__
#define __wuppertal_source_quark_smearing_h__

#include "meas/smear/quark_smearing.h"

namespace QDP
{
	extern Set TimeSliceSetQCDSF;
	void initTimeSliceSet();
}

namespace Chroma
{
  //! Name and registration
  /*! @ingroup smear */
  namespace WuppertalSourceQuarkSmearingEnv
  {
    extern const std::string name;
    bool registerAll ();
	void registerQDPGlobals();


    //! Params for Restricted Wuppertal quark smearing
    /*! @ingroup smear */
    struct Params
    {
      Params () {}
      Params (XMLReader& in, const std::string& path);
      void writeXML (XMLWriter& in, const std::string& path) const;

      Real wvf_param;                   /*!< Smearing width */
      int  wvfIntPar;                   /*!< Number of smearing hits */
      int  no_smear_dir;		/*!< No smearing in this direction */
      int  t_source;				/*!< only time slice that is smeared */
      bool smearSingleTimeslice;
      bool smearSingleSpin;
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
      QuarkSmear (const Params& p) : params (p) {
      	registerQDPGlobals();
      	tslice = params.smearSingleTimeslice ? QDP::TimeSliceSetQCDSF[params.t_source] : QDP::all ;
      }

      //! Smear the quark
      void operator() (T& quark, const multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      QuarkSmear () {}

    private:
      Params  params;   /*!< smearing params */
      Subset tslice;
    };

  }  // end namespace

  //! Reader
  /*! @ingroup smear */
  void read (XMLReader& xml, const std::string& path, WuppertalSourceQuarkSmearingEnv::Params& param);

  //! Writer
  /*! @ingroup smear */
  void write (XMLWriter& xml, const std::string& path, const WuppertalSourceQuarkSmearingEnv::Params& param);

}  // end namespace Chroma

#endif
