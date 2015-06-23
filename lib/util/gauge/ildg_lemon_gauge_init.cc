// $Id: ildg_lemon_gauge_init.cc,v 0.9 2011-05-16 22:52:41 bglaessle Exp $
/*! \file
 *  \brief Read a ILDG config with lemon
 */

#include "util/gauge/gauge_init_factory.h"
#include "util/gauge/gauge_init_aggregate.h"

#include "util/gauge/ildg_lemon_gauge_init.h"
#include "io/gauge_lemon_io.h"

#ifdef QDP_USE_LEMON

namespace Chroma
{

  // Read parameters
  void read(XMLReader& xml, const string& path, ILDGLemonGaugeInitEnv::Params& param)
  {
    ILDGLemonGaugeInitEnv::Params tmp(xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write(XMLWriter& xml, const string& path, const ILDGLemonGaugeInitEnv::Params& param)
  {
    param.writeXML(xml, path);
  }


  //! Hooks to register the class
  namespace ILDGLemonGaugeInitEnv
  {
    //! Callback function
    GaugeInit* createSource(XMLReader& xml_in,
			    const std::string& path)
    {
      return new GaugeIniter(Params(xml_in, path));
    }

    //! Name to be used
    const std::string name = "ILDG_LEMON";

    //! Local registration flag
    static bool registered = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::TheGaugeInitFactory::Instance().registerObject(name, createSource);
	registered = true;
      }
      return success;
    }


    // Parameters for running code
    Params::Params(XMLReader& xml, const string& path)
    {
      XMLReader paramtop(xml, path);

      read(paramtop, "cfg_file", cfg_file);
      // Default
      cfg_pario = QDPIO_PARALLEL;

      bool pario;
      if ( paramtop.count("ParallelIO") > 0 ) { 
        read(paramtop, "ParallelIO", pario);
        if( not pario ) { 
	  cfg_pario = QDPIO_SERIAL;
        }
      }
    }

    //! Parameters for running code
    void Params::writeXML(XMLWriter& xml, const string& path) const
    {
      push(xml, path);
    
      int version = 1;
      write(xml, "cfg_type", ILDGLemonGaugeInitEnv::name);
      write(xml, "cfg_file", cfg_file);
      if ( cfg_pario == QDPIO_PARALLEL ) { 
	bool pario = true;
	write(xml, "ParallelIO", pario);
      }
      pop(xml);
    }


    //! Returns a link smearing group with these params
    GroupXML_t   createXMLGroup(const Params& p)
    {
      GroupXML_t foo;

      XMLBufferWriter xml_tmp;
      write(xml_tmp, "Cfg", p);
      foo.xml = xml_tmp.printCurrentContext();
      foo.id = name;
      foo.path = "/Cfg";

      return foo;
    }


    // Initialize the gauge field
    void
    GaugeIniter::operator()(XMLReader& gauge_file_xml,
			    XMLReader& gauge_xml,
			    multi1d<LatticeColorMatrix>& u) const
    {
      u.resize(Nd);
      if( params.cfg_pario == QDPIO_SERIAL ) { 
	QDP_error_exit("Serial IO read on ILDG_LEMON not supported");
      }
      readGaugeILDGLemon(gauge_file_xml, gauge_xml, u, params.cfg_file, params.cfg_pario);
    }
  }
}

#endif // QDP_USE_LEMON

