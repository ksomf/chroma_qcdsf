/*! \file
 * \brief Inline calculation of the RMS of a wavefunction
 *
 * Calculates the root mean square of the wavefunction (from a propagator source)
 */

#include "meas/inline/hadron/inline_rms_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/inline/make_xml_file.h"
#include "meas/inline/io/named_objmap.h"
#include "meas/hadron/rms_w.h"
#include "io/qprop_io.h"

namespace Chroma
{
  namespace InlineRMSEnv
  {
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in,
					      const std::string& path)
      {
	return new InlineRMS(InlineRMSParams(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "RMS";

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
  } // end namespace


  //! RMS input
  void read(XMLReader& xml, const std::string& path, InlineRMSParams::NamedObject_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "source_id", input.source_id);
  }
  //! RMS input
  void read(XMLReader& xml, const std::string& path, InlineRMSParams::Param_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "version", input.version);
    read(inputtop, "psi_wfn", input.psi_wfn);
    read(inputtop, "psi_dag_psi_wfn", input.psi_dag_psi_wfn);
  }

  //! RMS output
  void write(XMLWriter& xml, const std::string& path, const InlineRMSParams::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "source_id", input.source_id);

    pop(xml);
  }
  //! RMS output
  void write(XMLWriter& xml, const std::string& path, const InlineRMSParams::Param_t& input)
  {
    push(xml, path);

    write(xml, "version", input.version);
    write(xml, "psi_wfn", input.psi_wfn);
    write(xml, "psi_dag_psi_wfn", input.psi_dag_psi_wfn);

    pop(xml);
  }


  // Param stuff
  InlineRMSParams::InlineRMSParams() { frequency = 0; }

  InlineRMSParams::InlineRMSParams(XMLReader& xml_in, const std::string& path)
  {
    try
    {
      XMLReader paramtop(xml_in, path);

      if (paramtop.count("Frequency") == 1)
	read(paramtop, "Frequency", frequency);
      else
	frequency = 1;

      read(paramtop, "Param", param);

      // Read in the source configuration info
      read(paramtop, "NamedObject", named_obj);

      // Possible alternate XML file pattern
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
  InlineRMSParams::writeXML(XMLWriter& xml_out, const std::string& path)
  {
    push(xml_out, path);

    write(xml_out, "Param", param);
    write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }


  // Function call
  void
  InlineRMS::operator()(unsigned long update_no,
			       XMLWriter& xml_out)
  {
    // If xml file not empty, then use alternate
    if (params.xml_file != "")
    {
		std::string xml_file = makeXMLFileName(params.xml_file, update_no);

      push(xml_out, "rms");
      write(xml_out, "update_no", update_no);
      write(xml_out, "xml_file", xml_file);
      pop(xml_out);

      XMLFileWriter xml(xml_file);
      func(update_no, xml);
    }
    else
    {
      func(update_no, xml_out);
    }
  }


  // Real work done here
  void
  InlineRMS::func(unsigned long update_no,
			 XMLWriter& xml_out)
  {
    START_CODE();

    QDPIO::cout << InlineRMSEnv::name << ": calculation of the RMS of a wavefunction" << std::endl;

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
      QDPIO::cerr << InlineRMSEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineRMSEnv::name << ": map call failed: " << e
		  << std::endl;
      QDP_abort(1);
    }
    const multi1d<LatticeColorMatrix>& u =
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "rms");
    write(xml_out, "update_no", update_no);

    proginfo(xml_out);    // Print out basic program info

    // Write out the input
    params.writeXML(xml_out, "Input");

    // Write out the config header
    write(xml_out, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 1);
    pop(xml_out);

    // Calculate some gauge invariant observables just for info.
    MesPlq(xml_out, "Observables", u);

    //
    // Read in the source along with relevant information.
    //
    XMLReader source_file_xml, source_record_xml;

    // Record the type of header
    bool make_sourceP = false;
    int j_decay;
    multi1d<int> srcloc(Nd);


    QDPIO::cout << "Snarf the source from a named buffer" << std::endl;
    try
    {
      // Try the cast to see if this is a valid source
      LatticePropagator& source_tmp =
	TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);

      // Snarf the source info. This is will throw if the source_id is not there
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getFileXML(source_file_xml);
      TheNamedObjMap::Instance().get(params.named_obj.source_id).getRecordXML(source_record_xml);

      // Try to invert this record XML into a source struct
      // First identify what kind of source might be here
      if (source_record_xml.count("/MakeSource") != 0)
      {
	make_sourceP = true;
	MakeSourceProp_t  orig_header;
	read(source_record_xml, "/MakeSource", orig_header);
	j_decay = orig_header.source_header.j_decay;
	srcloc = orig_header.source_header.getTSrce();
      }
      else
      {
	throw std::string("No appropriate header found");
      }

      // Write out the source header
      write(xml_out, "Source_file_info", source_file_xml);
      write(xml_out, "Source_record_info", source_record_xml);
    }
    catch (std::bad_cast)
    {
      QDPIO::cerr << InlineRMSEnv::name << ": caught dynamic cast error"
		  << std::endl;
      QDP_abort(1);
    }
    catch (const std::string& e)
    {
      QDPIO::cerr << InlineRMSEnv::name << ": error extracting source_header: " << e << std::endl;
      QDP_abort(1);
    }

    // Should be a valid cast now
    const LatticePropagator& quark_prop_source =
      TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.source_id);

    QDPIO::cout << "Source successfully read and parsed" << std::endl;

    // Sanity check - write out the norm2 of the source in the Nd-1 direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, Nd-1);

      multi1d<Double> source_corr = sumMulti(localNorm2(quark_prop_source),
					     phases.getSet());

      push(xml_out, "Source_correlator");
      write(xml_out, "source_corr", source_corr);
      pop(xml_out);
    }

    rms(quark_prop_source,j_decay,srcloc,params.param.psi_wfn,
	params.param.psi_dag_psi_wfn,xml_out,"rms_data");

    snoop.stop();
    QDPIO::cout << InlineRMSEnv::name << ": total time = "
		<< snoop.getTimeInSeconds()
		<< " secs" << std::endl;

    QDPIO::cout << InlineRMSEnv::name << ": ran successfully" << std::endl;

    END_CODE();
  }

}
