// $Id: inline_bda_w.cc, v 1.0 2011-04-11 chagen, rwschiel $

#include "meas/inline/hadron/inline_bda_w.h"
#include "meas/hadron/bdacomp_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/qprop_io.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"

namespace Chroma
{
  namespace InlineBDAEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement (XMLReader& xml_in,
					       const std::string& path)
      {
	return new InlineBDA (InlineBDAParams (xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "BDA";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance ().registerObject (name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! Reader for parameters
  void read (XMLReader& xml, const std::string& path, InlineBDAParams::Param_t& param)
  {
    XMLReader paramtop (xml, path);

    read (paramtop, "version", param.version);

    switch (param.version)
    {
    case 2:
      read (paramtop, "mom2_max", param.mom2_max);
      read (paramtop, "avg_equiv_mom", param.avg_equiv_mom);
      break;
    case 3: {
      multi1d<multi1d<int> > momdummy;
      read (paramtop, "moms", momdummy);
      param.moms.resize (momdummy.size (), 3);
      for (int i = 0; i < param.moms.size2 (); i++)
        param.moms [i] = momdummy [i];
      }
      break;
    default:
      QDPIO::cerr << "Input parameter version " << param.version << " unsupported." << std::endl;
      QDP_abort(1);
    }
    read (paramtop, "gam_diquark", param.gam_diquark);
  }

  //! Writer for parameters
  void write (XMLWriter& xml, const std::string& path, const InlineBDAParams::Param_t& param)
  {
    push (xml, path);

    write (xml, "version", param.version);

    switch (param.version) {
      case 2:
        write (xml, "mom2_max", param.mom2_max);
        write (xml, "avg_equiv_mom", param.avg_equiv_mom);
        break;
      case 3:
        multi1d<multi1d<int> > momdummy (param.moms.size2 ());
        for (int i = 0; i < momdummy.size (); i++) {
          momdummy [i].resize (param.moms.size1 ());
          momdummy [i] = param.moms [i];
        }
        write (xml, "moms", momdummy);
        break;
    }
    write (xml, "gam_diquark", param.gam_diquark);

    pop (xml);
  }



  //! Propagator input
  void read (XMLReader& xml, const std::string& path, InlineBDAParams::NamedObject_t& input)
  {
    XMLReader inputtop (xml, path);

    read (inputtop, "gauge_id", input.gauge_id);
    read (inputtop, "prop_ids", input.prop_ids);
  }

  //! Propagator output
  void write (XMLWriter& xml, const std::string& path, const InlineBDAParams::NamedObject_t& input)
  {
    push (xml, path);

    write (xml, "gauge_id", input.gauge_id);
    write (xml, "prop_ids", input.prop_ids);

    pop (xml);
  }



  // Param stuff
  InlineBDAParams::InlineBDAParams () {frequency = 0;}

  InlineBDAParams::InlineBDAParams (XMLReader& xml_in, const std::string& path)
  {
    try
    {
        XMLReader paramtop (xml_in, path);

        if (paramtop.count ("Frequency") == 1)
          read (paramtop, "Frequency", frequency);
        else
         frequency = 1;

        // Parameters for source construction
        read (paramtop, "Param", param);

        // Read in the output propagator/source configuration info
        read (paramtop, "NamedObject", named_obj);

        // Possible alternate qqq output file
        if (paramtop.count ("bda_file") != 0)
          read (paramtop, "bda_file", bda_file);

        // Possible alternate XML file pattern
        if (paramtop.count ("xml_file") != 0)
          read (paramtop, "xml_file", xml_file);

    }
    catch (const std::string& e)
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort (1);
    }
  }


  void
  InlineBDAParams::write (XMLWriter& xml_out, const std::string& path)
  {
    push (xml_out, path);

    Chroma::write (xml_out, "Param", param);
    Chroma::write (xml_out, "NamedObject", named_obj);
    QDP::write (xml_out, "bda_file", bda_file);
    QDP::write (xml_out, "xml_file", xml_file);

    pop (xml_out);
  }



  // Function call
  void
  InlineBDA::operator() (unsigned long update_no,
                         XMLWriter& xml_out)
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != ""){
		std::string xml_file = makeXMLFileName (params.xml_file, update_no);

      push (xml_out, "BDA_w");
      write (xml_out, "update_no", update_no);
      write (xml_out, "xml_file", xml_file);
      pop (xml_out);

      XMLFileWriter xml (xml_file);
      func (update_no, xml);
    }
    else
      func (update_no, xml_out);
  }


  // Real work done here
  void InlineBDA::func(unsigned long update_no,
		       XMLWriter& xml_out)
  {
    START_CODE ();

    StopWatch snoop;
    snoop.reset ();
    snoop.start ();


    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance ().getData<multi1d<LatticeColorMatrix> > (params.named_obj.gauge_id);
      TheNamedObjMap::Instance ().get (params.named_obj.gauge_id).getRecordXML (gauge_xml);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineBDAEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort (1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineBDAEnv::name << ": map call failed: " << e
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance ().getData<multi1d<LatticeColorMatrix> > (params.named_obj.gauge_id);

    push (xml_out, "BDA");
    write (xml_out, "update_no", update_no);
    QDPIO::cout << " BDA: Baryon distribution amplitudes for Wilson fermions" << std::endl;

    // Write out the input
    params.write (xml_out, "Input");

    proginfo (xml_out);    // Print out basic program info

    // Write out the input
    params.write (xml_out, "Input");


    // Write out the config info
    write (xml_out, "Config_info", gauge_xml);

    push (xml_out, "Output_version");
    write (xml_out, "out_version", 1);
    pop (xml_out);

    // First calculate some gauge invariant observables just for info.
    MesPlq (xml_out, "Observables", u);

    multi1d<ForwardProp_t> quark_header (params.named_obj.prop_ids.size ());
    multi1d<LatticePropagator> qprop (params.named_obj.prop_ids.size ());
    multi1d<Real> Mass (params.named_obj.prop_ids.size ());
    multi1d<std::string> sink_types (params.named_obj.prop_ids.size ());
    multi2d<int> bc (params.named_obj.prop_ids.size (), 4);

    // Now read the propagators we need
    for (int loop = 0; loop < params.named_obj.prop_ids.size (); loop++)
    {
      // Snarf the data into a copy
      qprop [loop] =
	TheNamedObjMap::Instance ().getData<LatticePropagator> (params.named_obj.prop_ids [loop]);

      // Snarf the source info.
      // This is will throw if the source_id is not there
      XMLReader prop_file_xml, prop_xml ;
      TheNamedObjMap::Instance ().get (params.named_obj.prop_ids [loop]).getFileXML (prop_file_xml);
      TheNamedObjMap::Instance ().get (params.named_obj.prop_ids [loop]).getRecordXML (prop_xml);
      // Try to invert this record XML into a ChromaProp struct
      // Also pull out the id of this source
      try
      {
	read (prop_xml, "/SinkSmear", quark_header [loop]);
	read (prop_xml, "/SinkSmear/PropSink/Sink/SinkType", sink_types [loop]);
      }
      catch (const std::string& e)
      {
	QDPIO::cerr << "Error extracting forward_prop header: " << e << std::endl;
	QDP_abort (1);
      }

      QDPIO::cout << "Try action and mass" << std::endl;
      Mass [loop] = getMass (quark_header [loop].prop_header.fermact);
      bc [loop]   = getFermActBoundary (quark_header [loop].prop_header.fermact);

      QDPIO::cout << "FermAct = " << quark_header [loop].prop_header.fermact.path << std::endl;
      QDPIO::cout << "Mass = " << Mass [loop] << std::endl;
      QDPIO::cout << "boundary = "
		  << bc [loop][0] << " "
		  << bc [loop][1] << " "
		  << bc [loop][2] << " "
		  << bc [loop][3] << std::endl;
    }


    // Derived from input prop
    int j_decay = quark_header [0].source_header.j_decay;
    multi1d<int> t_srce = quark_header [0].source_header.getTSrce ();
    int t_source = quark_header [0].source_header.t_source;
    int bc_spec = bc [0][j_decay] ;
    for (int loop = 0; loop < params.named_obj.prop_ids.size (); loop++)
    {
      if (quark_header [loop].source_header.j_decay != j_decay){
	QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << std::endl;
	QDP_abort(1);
      }
      if (bc [loop][j_decay] != bc_spec){
	QDPIO::cerr << "Error!! bc must be the same for all propagators " << std::endl;
	QDP_abort(1);
      }
      for (int d = 0; d < Nd; d++)
	if (quark_header [loop].source_header.t_source != t_source){
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << std::endl;
	  QDP_abort (1);
	}
    }

    // Initialize the slow Fourier transform phases
    SftMom* phases;
    switch (params.param.version) {
      case 2:
        phases = new SftMom (params.param.mom2_max, t_srce, params.param.avg_equiv_mom, j_decay);
        break;
      case 3:
        phases = new SftMom (params.param.moms, t_srce, j_decay);
        break;
    }

    // Write everything in binary. Maybe at some point one should switch to QDP.
    BinaryFileWriter bda_out (params.bda_file);

    multi2d<BaryonDA> qqq (phases->numMom (), phases->numSubsets ());

    // This routine computes the BaryonDAs.
    bdacomp (qqq, qprop [0], qprop [1], qprop [2], *phases, t_source, bc_spec, params.param.gam_diquark);

    write_bda (bda_out, qqq, *phases);

    bda_out.close ();

    snoop.stop ();
    QDPIO::cout << InlineBDAEnv::name << ": total time = "
		<< snoop.getTimeInSeconds ()
		<< " secs" << std::endl;

    QDPIO::cout << InlineBDAEnv::name << ": ran successfully" << std::endl;

    END_CODE ();
  }

}
