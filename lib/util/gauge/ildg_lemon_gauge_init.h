// -*- C++ -*-
// $Id: ildg_lemon_gauge_init.h,v 0.9 2012-05-16 22:06:42 bglaessle Exp $
/*! \file
 *  \brief Read a ILDG config with Lemon (mpi-IO)
 */

#ifndef __ildg_lemon_gauge_init_h__
#define __ildg_lemon_gauge_init_h__

#include "util/gauge/gauge_init.h"

#ifdef QDP_USE_LEMON

namespace Chroma
{

  //! Name and registration
  namespace ILDGLemonGaugeInitEnv
  {
    extern const std::string name;
    bool registerAll();
  

    //! Params for initializing config
    /*! @ingroup gauge */
    struct Params
    {
      Params() {}
      Params(XMLReader& in, const std::string& path);
      void writeXML(XMLWriter& in, const std::string& path) const;
    
      string cfg_file;		/*!< File name */
      QDP_serialparallel_t    cfg_pario;  /*!< QIO Parallel IO flag */
    };


    //! Returns a link smearing group with these params
    GroupXML_t   createXMLGroup(const Params& p);


    //! Gauge initialization
    /*! @ingroup gauge
     *
     * SZINQIO reader
     */
    class GaugeIniter : public GaugeInit
    {
    public:
      //! Full constructor
      GaugeIniter(const Params& p) : params(p) {}

      //! Initialize the gauge field
      void operator()(XMLReader& gauge_file_xml,
		      XMLReader& gauge_xml,
		      multi1d<LatticeColorMatrix>& u) const;

    private:
      //! Hide partial constructor
      GaugeIniter() {}

    private:
      Params  params;
    };

  }  // end namespace


  //! Reader
  /*! @ingroup gauge */
  void read(XMLReader& xml, const string& path, ILDGLemonGaugeInitEnv::Params& param);

  //! Writer
  /*! @ingroup gauge */
  void write(XMLWriter& xml, const string& path, const ILDGLemonGaugeInitEnv::Params& param);

}  // end namespace Chroma

#endif // QDP_USE_LEMON

#endif // __ildg_lemon_gauge_init_h__
