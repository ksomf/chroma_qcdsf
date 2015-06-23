// $Id: inline_hadspec_w.cc,v 3.16 2009-05-12 04:09:56 kostas Exp $
/*! \file
 * \brief Inline construction of hadron spectrum
 *
 * Spectrum calculations
 */

#include "meas/inline/hadron/inline_barspec_qcdsf_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "io/param_io.h"
#include "io/qprop_io.h"
#include "meas/hadron/barhqlq_qcdsf_w.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/smear/no_quark_displacement.h"

namespace Chroma
{
  namespace InlineBarSpecEnvQCDSF
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
					      const std::string& path)
      {
	return new InlineBarspecQCDSF(InlineBarspecParamsQCDSF(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "BARYON_SPECTRUM-QCDSF";

    //! Register all the factories
    bool registerAll()
    {
      bool success = true;
      if (! registered)
      {
	success &= TheInlineMeasurementFactory::Instance().registerObject(name, createMeasurement);
	registered = true;
      }
      return success;
    }
  }



  //! Reader for parameters
  void read(XMLReader& xml, const std::string& path, InlineBarspecParamsQCDSF::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version)
    {
    case 1:
      break;

    default:
      QDPIO::cerr << "Input parameter version " << version << " unsupported." << std::endl;
      QDP_abort(1);
    }


    if (paramtop.count("lime") != 0){
      read(paramtop, "lime", param.lime);
    } else
      param.lime = false;


    if (paramtop.count("xml") != 0){
      read(paramtop, "xml", param.xml);
    } else
      param.xml = true;



    read(paramtop, "time_rev", param.time_rev);
    if (paramtop.count("fwdbwd_average") != 0){
      read(paramtop, "fwdbwd_average", param.fwdbwd_average);
      param.fwdbwd_used = true;
    }
    else{
      param.fwdbwd_average = true;
      param.fwdbwd_used = false;
    }

    read(paramtop, "mom2_max", param.mom2_max);
    read(paramtop, "avg_equiv_mom", param.avg_equiv_mom);
  }


  //! Writer for parameters
  void write(XMLWriter& xml, const std::string& path, const InlineBarspecParamsQCDSF::Param_t& param)
  {
    push(xml, path);

    int version = 1;
    write(xml, "version", version);

    write(xml, "time_rev", param.time_rev);
    if(param.fwdbwd_used){
      write(xml, "fwdbwd_average", param.fwdbwd_average);
    }

    write(xml, "mom2_max", param.mom2_max);
    write(xml, "avg_equiv_mom", param.avg_equiv_mom);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineBarspecParamsQCDSF::NamedObject_t::Props_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "first_id", input.first_id);
    read(inputtop, "second_id", input.second_id);
    if (inputtop.count("third_id") != 0) {
      input.haveThird=true;
      read(inputtop, "third_id", input.third_id);
    } else {
      input.haveThird=false;
      input.third_id = input.second_id;
    }
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineBarspecParamsQCDSF::NamedObject_t::Props_t& input)
  {
    push(xml, path);

    write(xml, "first_id", input.first_id);
    write(xml, "second_id", input.second_id);
    if (input.second_id != input.third_id)
      write(xml, "third_id", input.third_id);

    pop(xml);
  }


  //! Propagator input
  void read(XMLReader& xml, const std::string& path, InlineBarspecParamsQCDSF::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "sink_pairs", input.sink_pairs);
  }

  //! Propagator output
  void write(XMLWriter& xml, const std::string& path, const InlineBarspecParamsQCDSF::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "sink_pairs", input.sink_pairs);

    pop(xml);
  }


  // Param stuff
  InlineBarspecParamsQCDSF::InlineBarspecParamsQCDSF()
  {
    frequency = 0;
  }

  InlineBarspecParamsQCDSF::InlineBarspecParamsQCDSF(XMLReader& xml_in, const std::string& path)
  {
    try
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      // Parameters for source construction
      read(paramtop, "Param", param);

      // Read in the output propagator/source configuration info
      read(paramtop, "NamedObject", named_obj);


      if (paramtop.count("lime_file") != 0)
      {
	read(paramtop, "lime_file", lime_file);
      }

      if (paramtop.count("xml_file") != 0)
      {
	read(paramtop, "xml_file", xml_file);
      }

    }
    catch(const std::string& e)
    {
      QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << std::endl;
      QDP_abort(1);
    }
  }


  void
  InlineBarspecParamsQCDSF::write(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);

    Chroma::write(xml_out, "Param", param);
    Chroma::write(xml_out, "NamedObject", named_obj);
    QDP::write(xml_out, "lime_file", lime_file);

    pop(xml_out);
  }



  // Anonymous namespace
  namespace
  {
    //! Useful structure holding sink props
    struct SinkPropContainer_t
    {
      ForwardProp_t prop_header;
	  std::string quark_propagator_id;
      Real Mass;

      multi1d<int> bc;

	  std::string source_type;
	  std::string source_disp_type;
	  std::string sink_type;
	  std::string sink_disp_type;
    };


    //! Useful structure holding sink props
    struct AllSinkProps_t
    {
      SinkPropContainer_t  sink_prop_1;
      SinkPropContainer_t  sink_prop_2;
      SinkPropContainer_t  sink_prop_3;
    };


    //! Read a sink prop
    void readSinkProp(SinkPropContainer_t& s, const std::string& id)
    {
      try
      {
	// Try a cast to see if it succeeds
	const LatticePropagator& foo =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(id);

	// Snarf the data into a copy
	s.quark_propagator_id = id;

	// Snarf the prop info. This is will throw if the prop_id is not there
	XMLReader prop_file_xml, prop_record_xml;
	TheNamedObjMap::Instance().get(id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(id).getRecordXML(prop_record_xml);

	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
		std::string xpath;
	  read(prop_record_xml, "/SinkSmear", s.prop_header);

	  read(prop_record_xml, "/SinkSmear/PropSource/Source/SourceType", s.source_type);
	  xpath = "/SinkSmear/PropSource/Source/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.source_disp_type);
	  else
	    s.source_disp_type = NoQuarkDisplacementEnv::getName();

	  read(prop_record_xml, "/SinkSmear/PropSink/Sink/SinkType", s.sink_type);
	  xpath = "/SinkSmear/PropSink/Sink/Displacement/DisplacementType";
	  if (prop_record_xml.count(xpath) != 0)
	    read(prop_record_xml, xpath, s.sink_disp_type);
	  else
	    s.sink_disp_type = NoQuarkDisplacementEnv::getName();
	}
      }
      catch( std::bad_cast )
      {
	QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": caught dynamic cast error"
		    << std::endl;
	QDP_abort(1);
      }
      catch (const std::string& e)
      {
	QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": error message: " << e
		    << std::endl;
	QDP_abort(1);
      }


      // Derived from input prop
      // Hunt around to find the mass
      // NOTE: this may be problematic in the future if actions are used with no
      // clear def. of a Mass
      QDPIO::cout << "Try action and mass" << std::endl;
      s.Mass = getMass(s.prop_header.prop_header.fermact);

      // Only baryons care about boundaries
      // Try to find them. If not present, assume dirichlet.
      // This turns off any attempt to time reverse which is the
      // only thing that the BC are affecting.
      s.bc.resize(Nd);
      s.bc = 0;

      try
      {
	s.bc = getFermActBoundary(s.prop_header.prop_header.fermact);
      }
      catch (const std::string& e)
      {
	QDPIO::cerr << InlineBarSpecEnvQCDSF::name
		    << ": caught exception. No BC found in these headers. Will assume dirichlet: " << e
		    << std::endl;
      }

      QDPIO::cout << "FermAct = " << s.prop_header.prop_header.fermact.id << std::endl;
      QDPIO::cout << "Mass = " << s.Mass << std::endl;
    }


    //! Read all sinks
    void readAllSinks(AllSinkProps_t& s,
		      InlineBarspecParamsQCDSF::NamedObject_t::Props_t sink_pair)
    {
      QDPIO::cout << "3rd id=" << sink_pair.third_id << std::endl;
      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.first_id << std::endl;
      readSinkProp(s.sink_prop_1, sink_pair.first_id);
      QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

      QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.second_id << std::endl;
      readSinkProp(s.sink_prop_2, sink_pair.second_id);
      QDPIO::cout << "Forward propagator successfully parsed" << std::endl;

      {
	QDPIO::cout << "Attempt to parse forward propagator = " << sink_pair.third_id << std::endl;
	readSinkProp(s.sink_prop_3, sink_pair.third_id);
	QDPIO::cout << "Forward propagator successfully parsed" << std::endl;
      }

    }

  } // namespace anonymous



  // Function call
  void
  InlineBarspecQCDSF::operator()(unsigned long update_no,
				 XMLWriter& xml_out)
  {
    if (params.param.lime)
      {
	if (params.lime_file != "")
	  {
	    push(xml_out, "barspec");
	    write(xml_out, "update_no", update_no);
	    write(xml_out, "lime_file", params.lime_file);
	    pop(xml_out);

	    func_lime(update_no, params.lime_file);
	  }
	else
	  {
	    QDPIO::cerr << "Error!! lime_file must be declared! " << std::endl;
	    QDP_abort(1);
	  }
      }

    if (params.param.xml)
      {
	if (params.xml_file != "")
	  {
		  std::string xml_file = makeXMLFileName(params.xml_file, update_no);

	    push(xml_out, "hadspec");
	    write(xml_out, "update_no", update_no);
	    write(xml_out, "xml_file", xml_file);
	    pop(xml_out);

	    XMLFileWriter xml(xml_file);
	    func_xml(update_no, xml);
	  }
	else
	  {
	    func_xml(update_no, xml_out);
	  }
      }


  }




  //
  // XML output
  //
  void
  InlineBarspecQCDSF::func_xml(unsigned long update_no,
			       XMLWriter& xml_out)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast )
    {
      QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": map call failed: " << e
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "hadspec");
    write(xml_out, "update_no", update_no);

    QDPIO::cout << " HADSPEC: Spectroscopy for Wilson-like fermions" << std::endl;
    QDPIO::cout << std::endl << "     Gauge group: SU(" << Nc << ")" << std::endl;
    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 15);
    pop(xml_out);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    // Keep an array of all the xml output buffers
    push(xml_out, "Wilson_hadron_measurements");

    // Now loop over the various fermion pairs
    for(int lpair=0; lpair < params.named_obj.sink_pairs.size(); ++lpair)
    {
      XMLBufferWriter xml_pair;

      push( xml_pair , "sink_pair" );

      const InlineBarspecParamsQCDSF::NamedObject_t::Props_t named_obj = params.named_obj.sink_pairs[lpair];

      push(xml_out, "elem");

      AllSinkProps_t all_sinks;
      readAllSinks(all_sinks, named_obj);

      // Derived from input prop
      multi1d<int> t_srce
                  = all_sinks.sink_prop_1.prop_header.source_header.getTSrce();
      int j_decay = all_sinks.sink_prop_1.prop_header.source_header.j_decay;
      int t0      = all_sinks.sink_prop_1.prop_header.source_header.t_source;

      // Sanity checks
      {
	if (all_sinks.sink_prop_2.prop_header.source_header.j_decay != j_decay)
	{
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_2.prop_header.source_header.t_source !=
	    all_sinks.sink_prop_1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.source_type != all_sinks.sink_prop_2.source_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.sink_type != all_sinks.sink_prop_2.sink_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
      }


      int bc_spec = 0;
      bc_spec = all_sinks.sink_prop_1.bc[j_decay] ;
      if (all_sinks.sink_prop_2.bc[j_decay] != bc_spec)
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}


      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, t_srce, params.param.avg_equiv_mom,
                    j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);

      // Masses
      write(xml_out, "Mass_1", all_sinks.sink_prop_1.Mass);
      write(xml_out, "Mass_2", all_sinks.sink_prop_2.Mass);
      if (named_obj.haveThird)
	write(xml_out, "Mass_3", all_sinks.sink_prop_3.Mass);
      write(xml_out, "t0", t0);

      // Save prop input
      push(xml_out, "Forward_prop_headers");
      write(xml_out, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
      write(xml_out, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
      if (named_obj.haveThird)
	write(xml_out, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
      pop(xml_out);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      push(xml_out, "Forward_prop_correlator");
      {
	const LatticePropagator& sink_prop_1 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);

	write(xml_out, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	write(xml_out, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	if (named_obj.haveThird)
	  write(xml_out, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));
      }
      pop(xml_out);


      push(xml_out, "SourceSinkType");
      {
	QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << std::endl;
	QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << std::endl;

	write(xml_out, "source_type_1", all_sinks.sink_prop_1.source_type);
	write(xml_out, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	write(xml_out, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	write(xml_out, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);

	write(xml_out, "source_type_2", all_sinks.sink_prop_2.source_type);
	write(xml_out, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	write(xml_out, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	write(xml_out, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);

	if (named_obj.haveThird) {
	  write(xml_out, "source_type_3", all_sinks.sink_prop_3.source_type);
	  write(xml_out, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	  write(xml_out, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	  write(xml_out, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	}
      }
      pop(xml_out);


      // References for use later
      const LatticePropagator& sink_prop_1 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
      const LatticePropagator& sink_prop_2 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
      const LatticePropagator& sink_prop_3 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);


      // Construct group name for output
	  std::string src_type;
      if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SF_POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "NORM_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "SF_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "RESMEAR_SOURCE-QCDSF")
 	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "SF_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "RAND_ZN_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "MOMENTUM_VOLUME_SOURCE")
        src_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDP_abort(1);
      }

	  std::string snk_type;
      if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	snk_type = "Point";
      else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "NORM_SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "RESMEAR_SINK-QCDSF")
	snk_type = "Shell";
     else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	snk_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDP_abort(1);
      }

	  std::string source_sink_type = src_type + "_" + snk_type;
      QDPIO::cout << "Source type = " << src_type << std::endl;
      QDPIO::cout << "Sink type = "   << snk_type << std::endl;

      barhqlq_qcdsf_xml(sink_prop_2, sink_prop_1, sink_prop_3, named_obj.haveThird, phases,
			t0, bc_spec, params.param.time_rev, params.param.fwdbwd_average,
			xml_out, source_sink_type + "_Wilson_Baryons");

      pop(xml_out);  // array element
    }
    pop(xml_out);  // Wilson_spectroscopy
    pop(xml_out);  // hadspec

    snoop.stop();
    QDPIO::cout << InlineBarSpecEnvQCDSF::name << ": total time = "
		<< snoop.getTimeInSeconds()
		<< " secs" << std::endl;

    QDPIO::cout << InlineBarSpecEnvQCDSF::name << ": ran successfully" << std::endl;

    END_CODE();
  }










  //
  // LIME output
  //
  void
  InlineBarspecQCDSF::func_lime(unsigned long update_no,
				std::string& lime_file)
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
    {
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
      TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
    }
    catch( std::bad_cast )
    {
      QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineBarSpecEnvQCDSF::name << ": map call failed: " << e
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    XMLBufferWriter xml_qcdsf;

    push(xml_qcdsf, "qcdsfDir");
    write(xml_qcdsf, "type", "barspecfn");
    write(xml_qcdsf, "update_no", update_no);

    QDPIO::cout <<         " BARSPEC-QCDSF: Spectroscopy for Wilson-like fermions" << std::endl;
    QDPIO::cout << std::endl << "            Gauge group: SU(" << Nc << ")" << std::endl;
    QDPIO::cout << "     volume: " << Layout::lattSize()[0];
    for (int i=1; i<Nd; ++i) {
      QDPIO::cout << " x " << Layout::lattSize()[i];
    }
    QDPIO::cout << std::endl;

    proginfo(xml_qcdsf);    // Print out basic program info

    // Write out the input
    params.write(xml_qcdsf, "Input");

    // Write out the config info
    write(xml_qcdsf, "Config_info", gauge_xml);

    push(xml_qcdsf, "Output_version");
    write(xml_qcdsf, "out_version", 15);
    pop(xml_qcdsf);


    // First calculate some gauge invariant observables just for info.
    MesPlq(xml_qcdsf, "Observables", u);

    // Keep an array of all the xml output buffers
    // push(xml_qcdsf, "Wilson_hadron_measurements");

    pop(xml_qcdsf);  // barspec



    QLimeWriter limewriter( params.lime_file.c_str() );

    QDPIO::cout << "writing LIME QCDSF header" << std::endl;
    uint64_t hdrsize = xml_qcdsf.str().length();
    limewriter.setRecordHeader( "qcdsfDir" , hdrsize , 1 , 0 );
    limewriter.write( (void *)( xml_qcdsf.str().c_str() ) , hdrsize );
    limewriter.endRecord();




    // Now loop over the various fermion pairs
    for(int lpair=0; lpair < params.named_obj.sink_pairs.size(); ++lpair)
    {
      bool lastSinkPair = (lpair == params.named_obj.sink_pairs.size()-1);
      XMLBufferWriter xml_pair;

      push( xml_pair , "sink_pair" );

      const InlineBarspecParamsQCDSF::NamedObject_t::Props_t named_obj = params.named_obj.sink_pairs[lpair];

      AllSinkProps_t all_sinks;
      readAllSinks(all_sinks, named_obj);

      // Derived from input prop
      multi1d<int> t_srce
                  = all_sinks.sink_prop_1.prop_header.source_header.getTSrce();
      int j_decay = all_sinks.sink_prop_1.prop_header.source_header.j_decay;
      int t0      = all_sinks.sink_prop_1.prop_header.source_header.t_source;

      // Sanity checks
      {
	if (all_sinks.sink_prop_2.prop_header.source_header.j_decay != j_decay)
	{
	  QDPIO::cerr << "Error!! j_decay must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_2.prop_header.source_header.t_source !=
	    all_sinks.sink_prop_1.prop_header.source_header.t_source)
	{
	  QDPIO::cerr << "Error!! t_source must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.source_type != all_sinks.sink_prop_2.source_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
	if (all_sinks.sink_prop_1.sink_type != all_sinks.sink_prop_2.sink_type)
	{
	  QDPIO::cerr << "Error!! source_type must be the same in a pair " << std::endl;
	  QDP_abort(1);
	}
      }


      int bc_spec = 0;

      bc_spec = all_sinks.sink_prop_1.bc[j_decay] ;
      if (all_sinks.sink_prop_2.bc[j_decay] != bc_spec)
	{
	  QDPIO::cerr << "Error!! bc must be the same for all propagators " << std::endl;
	  QDP_abort(1);
	}



      // Initialize the slow Fourier transform phases
      SftMom phases(params.param.mom2_max, t_srce, params.param.avg_equiv_mom,
                    j_decay);

      // Keep a copy of the phases with NO momenta
      SftMom phases_nomom(0, true, j_decay);

      // Masses
      write(xml_pair, "Mass_1", all_sinks.sink_prop_1.Mass);
      write(xml_pair, "Mass_2", all_sinks.sink_prop_2.Mass);
      if (named_obj.haveThird)
	write(xml_pair, "Mass_3", all_sinks.sink_prop_3.Mass);

      write(xml_pair, "t0", t0);

      // Save prop input
      push(xml_pair, "Forward_prop_headers");
      write(xml_pair, "First_forward_prop", all_sinks.sink_prop_1.prop_header);
      write(xml_pair, "Second_forward_prop", all_sinks.sink_prop_2.prop_header);
      if (named_obj.haveThird)
	write(xml_pair, "Third_forward_prop", all_sinks.sink_prop_3.prop_header);
      pop(xml_pair);

      // Sanity check - write out the norm2 of the forward prop in the j_decay direction
      // Use this for any possible verification
      push(xml_pair, "Forward_prop_correlator");
      {
	const LatticePropagator& sink_prop_1 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
	const LatticePropagator& sink_prop_2 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
	const LatticePropagator& sink_prop_3 =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);

	write(xml_pair, "forward_prop_corr_1", sumMulti(localNorm2(sink_prop_1), phases.getSet()));
	write(xml_pair, "forward_prop_corr_2", sumMulti(localNorm2(sink_prop_2), phases.getSet()));
	if (named_obj.haveThird)
	  write(xml_pair, "forward_prop_corr_3", sumMulti(localNorm2(sink_prop_3), phases.getSet()));

      }
      pop(xml_pair);


      push(xml_pair, "SourceSinkType");
      {
	QDPIO::cout << "Source_type_1 = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDPIO::cout << "Sink_type_1 = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDPIO::cout << "Source_type_2 = " << all_sinks.sink_prop_2.source_type << std::endl;
	QDPIO::cout << "Sink_type_2 = " << all_sinks.sink_prop_2.sink_type << std::endl;

	write(xml_pair, "source_type_1", all_sinks.sink_prop_1.source_type);
	write(xml_pair, "source_disp_type_1", all_sinks.sink_prop_1.source_disp_type);
	write(xml_pair, "sink_type_1", all_sinks.sink_prop_1.sink_type);
	write(xml_pair, "sink_disp_type_1", all_sinks.sink_prop_1.sink_disp_type);

	write(xml_pair, "source_type_2", all_sinks.sink_prop_2.source_type);
	write(xml_pair, "source_disp_type_2", all_sinks.sink_prop_2.source_disp_type);
	write(xml_pair, "sink_type_2", all_sinks.sink_prop_2.sink_type);
	write(xml_pair, "sink_disp_type_2", all_sinks.sink_prop_2.sink_disp_type);

	if (named_obj.haveThird) {
	  write(xml_pair, "source_type_3", all_sinks.sink_prop_3.source_type);
	  write(xml_pair, "source_disp_type_3", all_sinks.sink_prop_3.source_disp_type);
	  write(xml_pair, "sink_type_3", all_sinks.sink_prop_3.sink_type);
	  write(xml_pair, "sink_disp_type_3", all_sinks.sink_prop_3.sink_disp_type);
	}
      }
      pop(xml_pair);


      // References for use later
      const LatticePropagator& sink_prop_1 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_1.quark_propagator_id);
      const LatticePropagator& sink_prop_2 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_2.quark_propagator_id);
      const LatticePropagator& sink_prop_3 =
	TheNamedObjMap::Instance().getData<LatticePropagator>(all_sinks.sink_prop_3.quark_propagator_id);


      // Construct group name for output
	  std::string src_type;
      if (all_sinks.sink_prop_1.source_type == "POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SF_POINT_SOURCE")
	src_type = "Point";
      else if (all_sinks.sink_prop_1.source_type == "SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "NORM_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "SF_SHELL_SOURCE")
	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "RESMEAR_SOURCE-QCDSF")
 	src_type = "Shell";
      else if (all_sinks.sink_prop_1.source_type == "WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "SF_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "RAND_ZN_WALL_SOURCE")
	src_type = "Wall";
      else if (all_sinks.sink_prop_1.source_type == "MOMENTUM_VOLUME_SOURCE")
        src_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported source type = " << all_sinks.sink_prop_1.source_type << std::endl;
	QDP_abort(1);
      }

	  std::string snk_type;
      if (all_sinks.sink_prop_1.sink_type == "POINT_SINK")
	snk_type = "Point";
      else if (all_sinks.sink_prop_1.sink_type == "SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "NORM_SHELL_SINK")
	snk_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "RESMEAR_SINK-QCDSF")
 	src_type = "Shell";
      else if (all_sinks.sink_prop_1.sink_type == "WALL_SINK")
	snk_type = "Wall";
      else
      {
	QDPIO::cerr << "Unsupported sink type = " << all_sinks.sink_prop_1.sink_type << std::endl;
	QDP_abort(1);
      }

	  std::string source_sink_type = src_type + "_" + snk_type;
      QDPIO::cout << "Source type = " << src_type << std::endl;
      QDPIO::cout << "Sink type = "   << snk_type << std::endl;

      BaryonsQCDSF_t baryons;
      BaryonsQCDSF_t baryons_trev;
      barhqlq_qcdsf_lime(sink_prop_2, sink_prop_1, sink_prop_3, named_obj.haveThird, phases, t0,
			 bc_spec, params.param.time_rev, params.param.fwdbwd_average, baryons , baryons_trev);

      write(xml_pair,"baryon-number",baryons.barnum.size());

      pop( xml_pair );

      hdrsize = xml_pair.str().length();
      QDPIO::cout << "writing LIME record xml" << std::endl;
      limewriter.setRecordHeader( "meta-xml" , hdrsize , 0 , 0 );
      limewriter.write( (void *)( xml_pair.str().c_str() ) , hdrsize );
      limewriter.endRecord();

      // bool write_trev = !params.param.fwdbwd_average; // orig
      bool write_trev = params.param.time_rev &&!params.param.fwdbwd_average;
      {
	BinaryBufferWriter bar_bin;
	write( bar_bin , baryons );
	hdrsize = bar_bin.str().length();
	QDPIO::cout << "writing LIME binary data" << std::endl;
	limewriter.setRecordHeader( "baryons-bin" , hdrsize , 0 , (lastSinkPair && (!write_trev)) ? 1:0  );
	limewriter.write( (void *)( bar_bin.str().c_str() ) , hdrsize );
	limewriter.endRecord();
      }

      //if (!params.param.fwdbwd_average) // orig
      if(params.param.time_rev && !params.param.fwdbwd_average)
      {
	BinaryBufferWriter bar_bin;
	write( bar_bin , baryons_trev );
	hdrsize = bar_bin.str().length();
	QDPIO::cout << "writing LIME binary data" << std::endl;
	limewriter.setRecordHeader( "baryons-trev-bin" , hdrsize , 0 , lastSinkPair ? 1:0  );
	limewriter.write( (void *)( bar_bin.str().c_str() ) , hdrsize );
	limewriter.endRecord();
      }


      //pop(xml_pair);  // array element
    }

    //pop(xml_out);  // Wilson_spectroscopy

    //

    snoop.stop();
    QDPIO::cout << InlineBarSpecEnvQCDSF::name << ": total time = "
		<< snoop.getTimeInSeconds()
		<< " secs" << std::endl;

    QDPIO::cout << InlineBarSpecEnvQCDSF::name << ": ran successfully" << std::endl;

    END_CODE();
  }

};
