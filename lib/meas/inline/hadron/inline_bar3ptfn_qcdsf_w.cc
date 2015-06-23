// $Id: inline_bar3ptfn_w.cc,v 3.6 2007/06/10 14:40:23 edwards Exp $
/*! \file
 * \brief Inline measurement of bar3ptfn
 *
 * Form-factor measurements
 */

#include "meas/inline/hadron/inline_bar3ptfn_qcdsf_w.h"
//#include "meas/inline/hadron/inline_bar3ptfn_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "io/qprop_io.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/proginfo.h"
#include "meas/hadron/formfac_qcdsf_w.h"

#include "meas/inline/io/named_objmap.h"

namespace Chroma 
{ 
  namespace InlineBar3ptfnEnvQCDSF         
  { 
    namespace
    {
      AbsInlineMeasurement* createMeasurement(XMLReader& xml_in, 
					      const std::string& path) 
      {
	return new InlineBar3ptfnQCDSF(InlineBar3ptfnParamsQCDSF(xml_in, path));
      }

      //! Local registration flag
      bool registered = false;
    }

    const std::string name = "BAR3PTFN-QCDSF";

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


  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineBar3ptfnParamsQCDSF::SeqProp_t& input)
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "seqprop_id", input.seqprop_id);
    read(inputtop, "gamma_insertion", input.gamma_insertion);
  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineBar3ptfnParamsQCDSF::SeqProp_t& input)
  {
    push(xml, path);

    write(xml, "seqprop_id", input.seqprop_id);
    write(xml, "gamma_insertion", input.gamma_insertion);

    pop(xml);
  }


  //! Propagator parameters
  void read(XMLReader& xml, const string& path, InlineBar3ptfnParamsQCDSF::NamedObject_t& input , InlineBar3ptfnParamsQCDSF::Param_t& param )
  {
    XMLReader inputtop(xml, path);

    read(inputtop, "gauge_id", input.gauge_id);
    read(inputtop, "prop_id", input.prop_id);
    read(inputtop, "seqprops", input.seqprops);
    read(inputtop, "bar3ptfn_file", input.bar3ptfn_file);

    if (param.deriv >= 1)
      read(inputtop, "bar3ptfn_1D_file", input.bar3ptfn_1D_file);
    if (param.deriv >= 2)
      read(inputtop, "bar3ptfn_2D_file", input.bar3ptfn_2D_file);

  }

  //! Propagator parameters
  void write(XMLWriter& xml, const string& path, const InlineBar3ptfnParamsQCDSF::NamedObject_t& input)
  {
    push(xml, path);

    write(xml, "gauge_id", input.gauge_id);
    write(xml, "prop_id", input.prop_id);
    write(xml, "seqprops", input.seqprops);
    write(xml, "bar3ptfn_file", input.bar3ptfn_file);

    pop(xml);
  }


  // Reader for input parameters
  void read(XMLReader& xml, const string& path, InlineBar3ptfnParamsQCDSF::Param_t& param)
  {
    XMLReader paramtop(xml, path);

    int version;
    read(paramtop, "version", version);

    switch (version) 
      {
      case 6:
	// Uggh, assume j_decay = Nd-1 here. This could come from source.
	param.j_decay = Nd-1;
	break;

      case 7:
	read(paramtop, "j_decay", param.j_decay);
	break;

      case 8:
	read(paramtop, "j_decay", param.j_decay);
	read(paramtop, "deriv", param.deriv);
	break;

      default :
	/**************************************************************************/

	QDPIO::cerr << "Input parameter version " << version << " unsupported." << endl;
	QDP_abort(1);
      }

    read(paramtop, "mom2_max", param.mom2_max);
  }


  // Reader for input parameters
  /*  void write(XMLWriter& xml, const string& path, const InlineBar3ptfnParamsQCDSF::Param_t& param)
  {
    push(xml, path);

    int version = 6;

    write(xml, "version", version);
    write(xml, "mom2_max", param.mom2_max);

    pop(xml);
    }*/

  void write(XMLWriter& xml, const string& path, const InlineBar3ptfnParamsQCDSF::Param_t& param)
  {
    push(xml, path);
    int version = 8;
    write(xml, "version", version);
    write(xml, "mom2_max",param.mom2_max);
    write(xml, "j_decay",param.j_decay);
    write(xml, "deriv", param.deriv);
    pop(xml);
  }



  // Param stuff
  InlineBar3ptfnParamsQCDSF::InlineBar3ptfnParamsQCDSF()
  { 
    frequency = 0; 
  }

  InlineBar3ptfnParamsQCDSF::InlineBar3ptfnParamsQCDSF(XMLReader& xml_in, const std::string& path) 
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
	read(paramtop, "NamedObject", named_obj, param);
      }
    catch(const std::string& e) 
      {
	QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
	QDP_abort(1);
      }
  }


  void
  InlineBar3ptfnParamsQCDSF::write(XMLWriter& xml_out, const std::string& path) 
  {
    push(xml_out, path);
    
    // Parameters for source construction
    Chroma::write(xml_out, "Param", param);

    // Write out the output propagator/source configuration info
    Chroma::write(xml_out, "NamedObject", named_obj);

    pop(xml_out);
  }



  //--------------------------------------------------------------
  struct Output_versionQCDSF_t
  {
    int out_version;
  };

  struct FormFac_sequential_sourceQCDSF_t
  {
    string            seqsrc_type;
    int               t_source;
    int               t_sink;
    multi1d<int>      sink_mom;
    int               gamma_insertion;
    FormFac_insertionsQCDSF_t  formFacs;
  };

  struct FormFac_Wilson_3Pt_fn_measurementsQCDSF_t
  {
    int  output_version;   // Unique id for each output version of the structures
    multi1d<FormFac_sequential_sourceQCDSF_t> seqsrc;
  };

  struct Bar3ptfnQCDSF_t
  {
    Output_versionQCDSF_t                      output_version;
    InlineBar3ptfnParamsQCDSF::Param_t         param;
    FormFac_Wilson_3Pt_fn_measurementsQCDSF_t  bar;
  };




  void write(BinaryWriter& bin, const FormFac_sequential_sourceQCDSF_t& src)
  {
    //write(bin, src.seqsrc_type);
    //write(bin, src.t_source);
    //write(bin, src.t_sink);
    //write(bin, src.sink_mom);
    //write(bin, src.gamma_insertion);
    write(bin, src.formFacs);
  }
  void write(BinaryWriter& bin, const FormFac_Wilson_3Pt_fn_measurementsQCDSF_t& had)
  {
    //write(bin, had.output_version);
    write( bin , had.seqsrc , had.seqsrc.size() );
  }
  void write(BinaryWriter& bin, const Bar3ptfnQCDSF_t& bar)
  {
    //write(bin, "output_version", bar.output_version);
    //write(bin, "param", bar.param);
    write(bin, bar.bar);
  }



  void write(XMLWriter& xml, const string& path, const FormFac_momentaQCDSF_t& mom)
  {
    int magic = 20301;
    push(xml, path);
    write(xml, "magic",magic);
    write(xml, "inser_mom",mom.inser_mom);
    //write(xml, mom.local_current);
    //write(xml, mom.nonlocal_current);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const FormFac_insertionQCDSF_t& mes)
  {
    push(xml, path);
    write(xml, "gamma_value",mes.gamma_value);
    write(xml, "momenta",mes.momenta);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const FormFac_insertionsQCDSF_t& form)
  {
    push(xml, path);
    write(xml, "output_version",form.output_version);
    write(xml, "formFac",form.formFac);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const Output_versionQCDSF_t& ver)
  {
    push(xml, path);
    write(xml, "out_version",ver.out_version);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const FormFac_sequential_sourceQCDSF_t& src)
  {
    push(xml, path);
    write(xml, "seqsrc_type",src.seqsrc_type);
    write(xml, "t_source",src.t_source);
    write(xml, "t_sink",src.t_sink);
    write(xml, "sink_mom",src.sink_mom);
    write(xml, "gamma_insertion",src.gamma_insertion);
    write(xml, "formFacs",src.formFacs);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const FormFac_Wilson_3Pt_fn_measurementsQCDSF_t& had)
  {
    push(xml, path);
    write(xml, "output_version", had.output_version);
    write(xml, "seqsrc",had.seqsrc);
    pop(xml);
  }
  void write(XMLWriter& xml, const string& path, const Bar3ptfnQCDSF_t& bar)
  {
    push(xml, path);
    write(xml, "output_version", bar.output_version);
    write(xml, "param", bar.param);
    write(xml, "par", bar.bar);
    pop(xml);
  }

  void assertFermBC( const multi1d<int>& bc, const GroupXML_t& fermact ) {
    std::istringstream  xml_s(fermact.xml);
    XMLReader  fermacttop(xml_s);
    XMLReader  top(fermacttop, fermact.path);

    string bc_type;
    read( top, "descendant::FermBC[1]", bc_type );

    if ( bc_type != "SIMPLE_FERMBC" ) {
      QDPIO::cerr << __func__ << ": Currently bar3ptfn_qcsdf only supports SIMPLE_FERMBC.\nThis error might be unnecessary, if it is called because of information of propagators used to construct the sequential source,\n anyway check what you or the code is doing!" << endl;
      QDP_abort(1);
    }

    multi1d<int> bc_2 = getFermActBoundary( fermact );
    for( int mu=0; mu<Nd; ++mu )
        if ( bc[mu] != bc_2[mu] ) {
            QDPIO::cerr << __func__ << ": Currently bar3ptfn_qcsdf only supports matching boundary conditions.\nThis error might be unnecessary, if it is called because of information of propagators used to construct the sequential source,\n anyway check what you or the code is doing!" << endl;
            QDP_abort(1);
        }
  }


  // Function call
  void 
  InlineBar3ptfnQCDSF::operator()(unsigned long update_no,
				  XMLWriter& xml_out) 
  {
    START_CODE();

    StopWatch snoop;
    snoop.reset();
    snoop.start();

    XMLBufferWriter xml_qcdsf;

    push(xml_qcdsf, "qcdsfDir");

    write(xml_qcdsf, "type", "bar3ptfn");


    // Test and grab a reference to the gauge field
    XMLBufferWriter gauge_xml;
    try
      {
	TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);
	TheNamedObjMap::Instance().get(params.named_obj.gauge_id).getRecordXML(gauge_xml);
      }
    catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": map call failed: " << e 
		    << endl;
	QDP_abort(1);
      }
    const multi1d<LatticeColorMatrix>& u = 
      TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >(params.named_obj.gauge_id);

    push(xml_out, "bar3ptfn");
    write(xml_out, "update_no", update_no);
    push(xml_qcdsf, "bar3ptfn");
    write(xml_qcdsf, "update_no", update_no);

    QDPIO::cout << InlineBar3ptfnEnvQCDSF::name << ": Baryon form factors for Wilson fermions" << endl;

    proginfo(xml_out);    // Print out basic program info
    proginfo(xml_qcdsf);    // Print out basic program info

    // Write out the input
    params.write(xml_out, "Input");
    params.write(xml_qcdsf, "Input");

    // Write out the config info
    write(xml_out, "Config_info", gauge_xml);
    write(xml_qcdsf, "Config_info", gauge_xml);

    push(xml_out, "Output_version");
    write(xml_out, "out_version", 11);
    pop(xml_out);
    push(xml_qcdsf, "Output_version");
    write(xml_qcdsf, "out_version", 11);
    pop(xml_qcdsf);

    // First calculate some gauge invariant observables just for info.
    // This is really cheap.
    MesPlq(xml_out, "Observables", u);
    MesPlq(xml_qcdsf, "Observables", u);

    //
    // Read the quark propagator and extract headers
    //
    XMLReader prop_file_xml, prop_record_xml;
    LatticePropagator quark_propagator;
    ChromaProp_t prop_header;
    PropSourceConst_t source_header;
    QDPIO::cout << "Attempt to parse forward propagator" << endl;
    try
      {
	// Snarf the forward prop
	quark_propagator =
	  TheNamedObjMap::Instance().getData<LatticePropagator>(params.named_obj.prop_id);
	
	// Snarf the source info. This is will throw if the source_id is not there
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).getFileXML(prop_file_xml);
	TheNamedObjMap::Instance().get(params.named_obj.prop_id).getRecordXML(prop_record_xml);
   
	// Try to invert this record XML into a ChromaProp struct
	// Also pull out the id of this source
	{
	  read(prop_record_xml, "/Propagator/ForwardProp", prop_header);
	  read(prop_record_xml, "/Propagator/PropSource", source_header);
	}

	// Save propagator input
	write(xml_out, "Propagator_file_info", prop_file_xml);
	write(xml_out, "Propagator_record_info", prop_record_xml);
	write(xml_qcdsf, "Propagator_file_info", prop_file_xml);
	write(xml_qcdsf, "Propagator_record_info", prop_record_xml);
      }
    catch( std::bad_cast ) 
      {
	QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": caught dynamic cast error" 
		    << endl;
	QDP_abort(1);
      }
    catch (const string& e) 
      {
	QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": error message: " << e 
		    << endl;
	QDP_abort(1);
      }
    QDPIO::cout << "Forward propagator successfully parsed" << endl;

    // Derived from input prop
    multi1d<int> t_srce = source_header.getTSrce();
    int j_decay  = params.param.j_decay;
    int t_source = source_header.t_source;
    multi1d<int> bc = getFermActBoundary(prop_header.fermact);
    assertFermBC( bc, prop_header.fermact );

    // Sanity check - write out the norm2 of the forward prop in the j_decay direction
    // Use this for any possible verification
    {
      // Initialize the slow Fourier transform phases
      SftMom phases(0, true, j_decay);

      multi1d<Double> forward_prop_corr = sumMulti(localNorm2(quark_propagator), 
						   phases.getSet());

      push(xml_out, "Forward_prop_correlator");
      write(xml_out, "forward_prop_corr", forward_prop_corr);
      pop(xml_out);
      push(xml_qcdsf, "Forward_prop_correlator");
      write(xml_qcdsf, "forward_prop_corr", forward_prop_corr);
      pop(xml_qcdsf);
    }


    //
    // Big nested structure that is image of entire file
    //
    Bar3ptfnQCDSF_t  bar3pt, bar3pt_1D, bar3pt_2D;
    bar3pt.output_version.out_version = 11;  // bump this up everytime something changes
    bar3pt.param = params.param; // copy entire structure

    bar3pt_1D.output_version.out_version = 12;  // bump this up everytime something changes
    bar3pt_1D.param = params.param; // copy entire structure

    bar3pt_2D.output_version.out_version = 13;  // bump this up everytime something changes
    bar3pt_2D.param = params.param; // copy entire structure

    push(xml_out, "Wilson_3Pt_fn_measurements");
    push(xml_qcdsf, "Wilson_3Pt_fn_measurements");

    // Big nested structure that is image of all form-factors
    //    FormFac_Wilson_3Pt_fn_measurements_t  formfacs;
    bar3pt.bar.output_version = 4;  // bump this up everytime something changes
    bar3pt.bar.seqsrc.resize(params.named_obj.seqprops.size());

    bar3pt_1D.bar.output_version = 5;  // bump this up everytime something changes
    bar3pt_1D.bar.seqsrc.resize(params.named_obj.seqprops.size());

    bar3pt_2D.bar.output_version = 6;  // bump this up everytime something changes
    bar3pt_2D.bar.seqsrc.resize(params.named_obj.seqprops.size());

    XMLArrayWriter  xml_seq_src(xml_out, params.named_obj.seqprops.size());
    XMLArrayWriter  xml_seq_src_qcdsf(xml_qcdsf, params.named_obj.seqprops.size());
    push(xml_seq_src, "Sequential_source");
    push(xml_seq_src_qcdsf, "Sequential_source");

    for (int seq_src_ctr = 0; seq_src_ctr < params.named_obj.seqprops.size(); ++seq_src_ctr)
      {
	push(xml_seq_src);
	write(xml_seq_src, "seq_src_ctr", seq_src_ctr);
	push(xml_seq_src_qcdsf);
	write(xml_seq_src_qcdsf, "seq_src_ctr", seq_src_ctr);

	// Read the sequential propagator
	// Read the quark propagator and extract headers
	LatticePropagator seq_quark_prop;
	SeqSource_t seqsource_header;
	QDPIO::cout << "Attempt to parse sequential propagator" << endl;
	try
	  {
	    std::string seqprop_id = params.named_obj.seqprops[seq_src_ctr].seqprop_id;

	    // Snarf the backward prop
	    seq_quark_prop =
	      TheNamedObjMap::Instance().getData<LatticePropagator>(seqprop_id);
	
	    // Snarf the source info. This is will throw if the source_id is not there
	    XMLReader seqprop_file_xml, seqprop_record_xml;
	    TheNamedObjMap::Instance().get(seqprop_id).getFileXML(seqprop_file_xml);
	    TheNamedObjMap::Instance().get(seqprop_id).getRecordXML(seqprop_record_xml);
   
	    // Try to invert this record XML into a ChromaProp struct
	    // Also pull out the id of this source
	    // NEED SECURITY HERE - need a way to cross check props. Use the ID.
	    {
	      read(seqprop_record_xml, "/SequentialProp/SeqSource", seqsource_header);


              // code to catch "bad" boundary conditions ...
              SequentialProp_t seq_header;
	      read(seqprop_record_xml, "/SequentialProp", seq_header);

              // - absolutely necessary!
              assertFermBC( bc, seq_header.seqprop_header.fermact );
              
              // - this might be physically ok, in the unlikely case that it happens
              //   the corresponding person should better check it again
              for( int j=0; j<seq_header.forward_props.size(); ++j )
                assertFermBC( bc, seq_header.forward_props[j].prop_header.fermact );
	    }

	    // Save seqprop input
	    write(xml_seq_src, "SequentialProp_file_info", seqprop_file_xml);
	    write(xml_seq_src, "SequentialProp_record_info", seqprop_record_xml);
	    write(xml_seq_src_qcdsf, "SequentialProp_file_info", seqprop_file_xml);
	    write(xml_seq_src_qcdsf, "SequentialProp_record_info", seqprop_record_xml);
	  }
	catch( std::bad_cast ) 
	  {
	    QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": caught dynamic cast error" 
			<< endl;
	    QDP_abort(1);
	  }
	catch (const string& e) 
	  {
	    QDPIO::cerr << InlineBar3ptfnEnvQCDSF::name << ": map call failed: " << e 
			<< endl;
	    QDP_abort(1);
	  }
	QDPIO::cout << "Sequential propagator successfully parsed" << endl;

	// Sanity check - write out the norm2 of the forward prop in the j_decay direction
	// Use this for any possible verification
	{
	  // Initialize the slow Fourier transform phases
	  SftMom phases(0, true, Nd-1);
      
	  multi1d<Double> backward_prop_corr = sumMulti(localNorm2(seq_quark_prop), 
							phases.getSet());
      
	  push(xml_seq_src, "Backward_prop_correlator");
	  write(xml_seq_src, "backward_prop_corr", backward_prop_corr);
	  pop(xml_seq_src);
	  push(xml_seq_src_qcdsf, "Backward_prop_correlator");
	  write(xml_seq_src_qcdsf, "backward_prop_corr", backward_prop_corr);
	  pop(xml_seq_src_qcdsf);
	}

	// Use extra gamma insertion
	int gamma_insertion = params.named_obj.seqprops[seq_src_ctr].gamma_insertion;

	// Derived from input seqprop
	std::string   seqsrc_type = seqsource_header.seqsrc.id;
	QDPIO::cout << "Seqsource name = " << seqsrc_type  << endl;
	int           t_sink   = seqsource_header.t_sink;
	multi1d<int>  sink_mom = seqsource_header.sink_mom;

	write(xml_seq_src, "hadron_type", "HADRON");
	write(xml_seq_src, "seqsrc_type", seqsrc_type);
	write(xml_seq_src, "t_source", t_source);
	write(xml_seq_src, "t_sink", t_sink);
	write(xml_seq_src, "sink_mom", sink_mom);
	write(xml_seq_src, "gamma_insertion", gamma_insertion);

	write(xml_seq_src_qcdsf, "hadron_type", "HADRON");
	write(xml_seq_src_qcdsf, "seqsrc_type", seqsrc_type);
	write(xml_seq_src_qcdsf, "t_source", t_source);
	write(xml_seq_src_qcdsf, "t_sink", t_sink);
	write(xml_seq_src_qcdsf, "sink_mom", sink_mom);
	write(xml_seq_src_qcdsf, "gamma_insertion", gamma_insertion);

	bar3pt.bar.seqsrc[seq_src_ctr].seqsrc_type   = seqsrc_type;
	bar3pt.bar.seqsrc[seq_src_ctr].t_source      = t_source;
	bar3pt.bar.seqsrc[seq_src_ctr].t_sink        = t_sink;
	bar3pt.bar.seqsrc[seq_src_ctr].sink_mom      = sink_mom;
	bar3pt.bar.seqsrc[seq_src_ctr].gamma_insertion = gamma_insertion;

	bar3pt_1D.bar.seqsrc[seq_src_ctr].seqsrc_type   = seqsrc_type;
	bar3pt_1D.bar.seqsrc[seq_src_ctr].t_source      = t_source;
	bar3pt_1D.bar.seqsrc[seq_src_ctr].t_sink        = t_sink;
	bar3pt_1D.bar.seqsrc[seq_src_ctr].sink_mom      = sink_mom;
	bar3pt_1D.bar.seqsrc[seq_src_ctr].gamma_insertion = gamma_insertion;

	bar3pt_2D.bar.seqsrc[seq_src_ctr].seqsrc_type   = seqsrc_type;
	bar3pt_2D.bar.seqsrc[seq_src_ctr].t_source      = t_source;
	bar3pt_2D.bar.seqsrc[seq_src_ctr].t_sink        = t_sink;
	bar3pt_2D.bar.seqsrc[seq_src_ctr].sink_mom      = sink_mom;
	bar3pt_2D.bar.seqsrc[seq_src_ctr].gamma_insertion = gamma_insertion;
	
        multi1d <LatticeColorMatrix> u_deriv=u; //! To correct the Nabla BUG
        for ( int dir=0; dir<Nd; ++dir ) {
            if( bc[dir]!=1 ) {
                QDPIO::cout << "Applying BC " << bc[dir] 
                            << " on gauge field_" << dir << " for Nabla"<< endl;
                u_deriv[dir] *= where( 
                        Layout::latticeCoordinate(dir)==Layout::lattSize()[dir]-1,
                        Double(bc[dir]), 
                        Double(1.) );
            }
        }

	// Now the 3pt contractions
	SftMom phases(params.param.mom2_max, t_srce, sink_mom, false, j_decay);

	FormFacQCDSF(bar3pt.bar.seqsrc[seq_src_ctr].formFacs, 
		     u_deriv, quark_propagator, seq_quark_prop, gamma_insertion,
		     phases, t_source);

	if (params.param.deriv >= 1) {
	  FormFac1DQCDSF(bar3pt_1D.bar.seqsrc[seq_src_ctr].formFacs, 
			 u_deriv, quark_propagator, seq_quark_prop, gamma_insertion,
			 phases, t_source);
	}

	if (params.param.deriv >= 2) {
	  FormFac2DQCDSF(bar3pt_2D.bar.seqsrc[seq_src_ctr].formFacs, 
			 u_deriv, quark_propagator, seq_quark_prop, gamma_insertion,
			 phases, t_source);
	}


	pop(xml_seq_src);   // elem
	pop(xml_seq_src_qcdsf);   // elem
      } // end loop over sequential sources

    pop(xml_seq_src);  // Sequential_source
    pop(xml_seq_src_qcdsf);  // Sequential_source

    //    BinaryWFileriter  bin_out(params.named_obj.bar3ptfn_file);
    //    write(bin_out, bar3ptfn);
    //    bin_out.close();

    pop(xml_out);  // Wilson_3Pt_fn_measurements
    pop(xml_qcdsf);  // Wilson_3Pt_fn_measurements

    // Close the namelist output file XMLDAT
    pop(xml_out);     // bar3ptfn
    pop(xml_qcdsf);     // bar3ptfn

    pop(xml_qcdsf);     // qcdsf

    // no derivative
    {
      QLimeWriter w( params.named_obj.bar3ptfn_file.c_str() );

      uint64_t hdrsize = xml_qcdsf.str().length();
      w.setRecordHeader( "meta-xml" , hdrsize , 1 , 0 );
      w.write( (void *)( xml_qcdsf.str().c_str() ) , hdrsize );
      w.endRecord();

      XMLBufferWriter bar_xml;
      push(bar_xml,"bar3ptfn");
      write( bar_xml , "records" , bar3pt );
      write( bar_xml , "deriv" , 0 );
      pop(bar_xml);
      hdrsize = bar_xml.str().length();
      w.setRecordHeader( "bar3ptfn-xml" , hdrsize , 0 , 0 );
      w.write( (void *)( bar_xml.str().c_str() ) , hdrsize );
      w.endRecord();

      BinaryBufferWriter bar_bin;
      write( bar_bin , bar3pt );
      hdrsize = bar_bin.str().length();
      w.setRecordHeader( "bar3ptfn-bin" , hdrsize , 0 , 1 );
      w.write( (void *)( bar_bin.str().c_str() ) , hdrsize );
      w.endRecord();
    }

    // 1D
    if (params.param.deriv >= 1)
      {
	QLimeWriter w( params.named_obj.bar3ptfn_1D_file.c_str() );

	uint64_t hdrsize = xml_qcdsf.str().length();
	w.setRecordHeader( "meta-xml" , hdrsize , 1 , 0 );
	w.write( (void *)( xml_qcdsf.str().c_str() ) , hdrsize );
	w.endRecord();

	XMLBufferWriter bar_xml;
	push(bar_xml,"bar3ptfn");
	write( bar_xml , "records" , bar3pt_1D );
	write( bar_xml , "deriv" , 1 );
	pop(bar_xml);
	hdrsize = bar_xml.str().length();
	w.setRecordHeader( "bar3ptfn-xml" , hdrsize , 0 , 0 );
	w.write( (void *)( bar_xml.str().c_str() ) , hdrsize );
	w.endRecord();

	BinaryBufferWriter bar_bin;
	write( bar_bin , bar3pt_1D );
	hdrsize = bar_bin.str().length();
	w.setRecordHeader( "bar3ptfn-bin" , hdrsize , 0 , 1 );
	w.write( (void *)( bar_bin.str().c_str() ) , hdrsize );
	w.endRecord();
      }

    // 2D
    if (params.param.deriv >= 2)
      {
	QLimeWriter w( params.named_obj.bar3ptfn_2D_file.c_str() );

	uint64_t hdrsize = xml_qcdsf.str().length();
	w.setRecordHeader( "meta-xml" , hdrsize , 1 , 0 );
	w.write( (void *)( xml_qcdsf.str().c_str() ) , hdrsize );
	w.endRecord();

	XMLBufferWriter bar_xml;
	push(bar_xml,"bar3ptfn");
	write( bar_xml , "records" , bar3pt_2D );
	write( bar_xml , "deriv" , 2 );
	pop(bar_xml);
	hdrsize = bar_xml.str().length();
	w.setRecordHeader( "bar3ptfn-xml" , hdrsize , 0 , 0 );
	w.write( (void *)( bar_xml.str().c_str() ) , hdrsize );
	w.endRecord();

	BinaryBufferWriter bar_bin;
	write( bar_bin , bar3pt_2D );
	hdrsize = bar_bin.str().length();
	w.setRecordHeader( "bar3ptfn-bin" , hdrsize , 0 , 1 );
	w.write( (void *)( bar_bin.str().c_str() ) , hdrsize );
	w.endRecord();
      }

    snoop.stop();
    QDPIO::cout << InlineBar3ptfnEnvQCDSF::name << ": total time = "
		<< snoop.getTimeInSeconds() 
		<< " secs" << endl;

    QDPIO::cout << InlineBar3ptfnEnvQCDSF::name << ": ran successfully" << endl;

    END_CODE();
  } 

};
