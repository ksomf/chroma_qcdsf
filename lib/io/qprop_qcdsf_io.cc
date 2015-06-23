// $Id: qprop_qcdsf_io.cc, v 1.0 2011-08-31 15:26:29 bglaessle Exp $
/*! \file
 * \brief Modified structs for source / sink Metadata, for resmeared objects
 */

//#define DEBUG_QPROP_QCDSF_IO

#include "chromabase.h"
#include "io/param_io.h"
#include "io/qprop_qcdsf_io.h"

#include "meas/smear/simple_quark_displacement.h"
#include "meas/smear/ape_link_smearing.h"

namespace Chroma { 
	
	void readSmearingArray(XMLReader& xml, const string& path, multi1d<GroupXML_t>& arr ) {
		XMLReader paramtop(xml, path);
		
		int num_smearings= paramtop.count( "elem" );
		arr.resize( num_smearings );

		for( int k=0; k<num_smearings; ++k ) {
			stringstream sstr;
			sstr << "elem[" << k+1 << "]/SmearingParam";
			arr[k] = readXMLGroup( paramtop, sstr.str(), "wvf_kind" );
		}
	}

	void writeSmearingArray(XMLWriter& xml, const string& path, const multi1d<GroupXML_t>& arr ) {
		push( xml, path );
		for( int k=0; k<arr.size(); ++k ) {
			push( xml, "elem" );
			xml << arr[k].xml;
			pop( xml ); // "elem"
		}
		pop( xml ); // path
	}

	// ResmearSourcePropQCDSF_t reader
	void read(XMLReader& xml, const string& path, ResmearSourcePropQCDSF_t& param) {
		XMLReader paramtop(xml, path);

		string source_type;
		read( paramtop, "PropSource/Source/SourceType", source_type );
		
		if( source_type=="RESMEAR_SOURCE-QCDSF" ) {
			QDPIO::cout << __func__ << " ResmearSourcePropQCDSF_t: resmearing RESMEAR_SOURCE-QCDSF" << endl;
			read(paramtop, "PropSource/Source/OriginalSource", param.orig_source_header );

			readSmearingArray(paramtop, "PropSource/Source/SmearingHeaders", param.smearing_headers );
			param.fake_smear_header = readXMLGroup(paramtop, "PropSource/Source/SmearingParam", "wvf_kind");
		}
		else {		
			QDPIO::cout << __func__ << " ResmearSourcePropQCDSF_t: resmearing MAKE_SOURCE" << endl;

		//	Reading in Previous MAKE_SOURCE object
			read(paramtop, "PropSource", param.orig_source_header );

			if ( paramtop.count("PropSource/Source/SmearingParam")==1 ) {
				param.fake_smear_header = readXMLGroup(paramtop, "PropSource/Source/SmearingParam", "wvf_kind");
			//	param.smearing_headers.resize(1);
			//	param.smearing_headers[0] = param.fake_smear_header; // use as initial
			}
			else
				QDPIO::cout << __func__ << " ResmearSourcePropQCDSF_t: No previous smearing found" << endl;
		}

		readGaugeHeader(paramtop, "Config_info", param.gauge_header);

#ifdef DEBUG_QPROP_QCDSF_IO
		QDPIO::cout << "XX" << param.orig_source_header.source.xml << "XX" << endl;
		QDPIO::cout << "XX" << param.fake_smear_header.xml << "XX" << endl;
		QDPIO::cout << "XX" << param.smearing_headers.size() << "XX" << endl;
		QDPIO::cout << "XX" << param.smearing_headers[0].xml << "XX" << endl;
		QDPIO::cout << "XX" << param.gauge_header << "XX" << endl;
#endif
	}

	// ResmearSourcePropQCDSF_t writer
	void write(XMLWriter& xml, const string& path, const ResmearSourcePropQCDSF_t& param) {
		push(xml, path);

		push(xml, "PropSource");
		
		int version = 6; // faking a proper Chroma XML structure version 6
		write(xml, "version", version);
		
		push(xml, "Source");

		version = 1; // RESMEAR Version
	//	string source_type = "SHELL_SOURCE";
		string source_type = "RESMEAR_SOURCE-QCDSF";
		write(xml, "version", version );
		write(xml, "SourceType", source_type );
		write(xml, "t_source", param.orig_source_header.t_source );
		write(xml, "j_decay", param.orig_source_header.j_decay );

		// additional elements, which are necessary for successive stuff
		write(xml, "OriginalSource", param.orig_source_header );
		writeSmearingArray(xml, "SmearingHeaders", param.smearing_headers);

		// FIX SMEARING
		{
			int wvfIntPar_1, wvfIntPar_2;
			Real wvf_param_1, wvf_param_2;
			int no_smear_dir_1, no_smear_dir_2;

			{
				int this_smear = param.smearing_headers.size()-1;
				std::istringstream  xml_s2(param.smearing_headers[this_smear].xml );
				XMLReader inputtop2(xml_s2);
				XMLReader smeartop2(inputtop2, "/SmearingParam" );

				read(smeartop2,"wvfIntPar", wvfIntPar_2);
				read(smeartop2,"wvf_param", wvf_param_2);
				read(smeartop2,"no_smear_dir", no_smear_dir_2);
			}

			if( param.fake_smear_header.xml==""
				or param.fake_smear_header.id==""
				or param.fake_smear_header.path=="" ) {
				wvfIntPar_1 = 0;
				wvf_param_1 = 0;
				no_smear_dir_1 = no_smear_dir_2;
			}
			else {
				std::istringstream  xml_s(param.fake_smear_header.xml );
				XMLReader inputtop(xml_s);
				XMLReader smeartop(inputtop, "/SmearingParam" );

				read(smeartop,"wvfIntPar", wvfIntPar_1 );
				read(smeartop,"wvf_param", wvf_param_1 );
				read(smeartop,"no_smear_dir", no_smear_dir_1 );
			}
			
			string wvf_kind = "MIXED_QUARK_SMEARING";
			int wvfIntPar = wvfIntPar_1 + wvfIntPar_2;
			Real wvf_param = ( wvfIntPar_1*wvf_param_1 + wvfIntPar_2*wvf_param_2 ) / Real(wvfIntPar);
			int no_smear_dir = ( no_smear_dir_1 == no_smear_dir_2 ? no_smear_dir_1 : -2 );

			push(xml, "SmearingParam" );
		
			write(xml, "wvf_kind", wvf_kind );
			write(xml, "wvf_param", wvf_param );
			write(xml, "wvfIntPar", wvfIntPar );
			write(xml, "no_smear_dir", no_smear_dir );
		
			pop(xml); // SmearingParam
		}

		pop(xml); // Source
		pop(xml); // PropSource

	//	writeModifiedSourceHeader(xml, "PropSource", param);
	//	write(xml, "PropSource", param.orig_source_header);
	//	write(xml, "SmearingHeaders", param.smearing_headers);
		write(xml, "Config_info", param.gauge_header);

		pop(xml);		
	}

	// ResmearForwardPropQCDSF_t reader
	void read(XMLReader& xml, const string& path, ResmearForwardPropQCDSF_t& param) {
		XMLReader paramtop(xml, path);

		string sink_type;
		read( paramtop, "PropSink/Sink/SinkType", sink_type );

		if( sink_type=="RESMEAR_SINK-QCDSF" ) {
			QDPIO::cout << __func__ << " ResmearForwardPropQCDSF_t: resmearing RESMEAR_SINK-QCDSF" << endl;
			read(paramtop, "PropSink/Sink/OriginalSink", param.orig_sink_header );

			readSmearingArray(paramtop, "PropSink/Sink/SmearingHeaders", param.smearing_headers );
			param.fake_smear_header = readXMLGroup(paramtop, "PropSink/Sink/SmearingParam", "wvf_kind");
		}
		else {		
			QDPIO::cout << __func__ << " ResmearForwardPropQCDSF_t: resmearing SINK_SMEAR" << endl;

		//	Reading in Previous SINK_SMEAR object
			read(paramtop, "PropSink", param.orig_sink_header );

			if ( paramtop.count("PropSink/Sink/SmearingParam")==1 ) {
				param.fake_smear_header = readXMLGroup(paramtop, "PropSink/Sink/SmearingParam", "wvf_kind");
			//	param.smearing_headers.resize(1);
			//	param.smearing_headers[0] = param.fake_smear_header; // use as initial
			}
			else
				QDPIO::cout << __func__ << " ResmearForwardPropQCDSF_t: No previous smearing found" << endl;
		}

		read(paramtop, "ForwardProp", param.prop_header);
		read(paramtop, "PropSource", param.source_header);
		readGaugeHeader(paramtop, "Config_info", param.gauge_header);
	}

	// ForwardPropQCDSFt writer
	void write(XMLWriter& xml, const string& path, const ResmearForwardPropQCDSF_t& param) {
		// if( path != "." )
		push(xml, path);

//		int version = 1;
//		write(xml, "version", version);

		push(xml, "PropSink");
		
		int version = 5; // faking a proper Chroma XML structure version ???
		write(xml, "version", version);
		
		push(xml, "Sink");

		version = 1; // RESMEAR Version
	//	string source_type = "SHELL_SINK";
		string source_type = "RESMEAR_SINK-QCDSF";
		write(xml, "version", version );
		write(xml, "SinkType", source_type );
		write(xml, "j_decay", param.orig_sink_header.j_decay );

		// additional elements, which are necessary for successive stuff
		write(xml, "OriginalSink", param.orig_sink_header );
		writeSmearingArray(xml, "SmearingHeaders", param.smearing_headers);

		{
			// FIX SMEARING
			int wvfIntPar_1, wvfIntPar_2;
			Real wvf_param_1, wvf_param_2;
			int no_smear_dir_1, no_smear_dir_2;

			{
				int this_smear = param.smearing_headers.size()-1;
				std::istringstream  xml_s2(param.smearing_headers[this_smear].xml );
				XMLReader inputtop2(xml_s2);
				XMLReader smeartop2(inputtop2, "/SmearingParam" );

				read(smeartop2,"wvfIntPar", wvfIntPar_2);
				read(smeartop2,"wvf_param", wvf_param_2);
				read(smeartop2,"no_smear_dir", no_smear_dir_2);
			}

			if( param.fake_smear_header.xml==""
				or param.fake_smear_header.id==""
				or param.fake_smear_header.path=="" ) {
				wvfIntPar_1 = 0;
				wvf_param_1 = 0;
				no_smear_dir_1 = no_smear_dir_2;
			}
			else {
				std::istringstream  xml_s(param.fake_smear_header.xml );
				XMLReader inputtop(xml_s);
				XMLReader smeartop(inputtop, "/SmearingParam" );

				read(smeartop,"wvfIntPar", wvfIntPar_1 );
				read(smeartop,"wvf_param", wvf_param_1 );
				read(smeartop,"no_smear_dir", no_smear_dir_1 );
			}
			
			string wvf_kind = "MIXED_QUARK_SMEARING";
			int wvfIntPar = wvfIntPar_1 + wvfIntPar_2;
			Real wvf_param = ( wvfIntPar_1*wvf_param_1 + wvfIntPar_2*wvf_param_2 ) / Real(wvfIntPar);
			int no_smear_dir = ( no_smear_dir_1 == no_smear_dir_2 ? no_smear_dir_1 : -2 );

			push(xml, "SmearingParam" );
		
			write(xml, "wvf_kind", wvf_kind );
			write(xml, "wvf_param", wvf_param );
			write(xml, "wvfIntPar", wvfIntPar );
			write(xml, "no_smear_dir", no_smear_dir );
		
			pop(xml); // SmearingParam
		}

		pop(xml); // Sink
		pop(xml); // PropSink

		write(xml, "ForwardProp", param.prop_header);
		write(xml, "PropSource", param.source_header);
		write(xml, "Config_info", param.gauge_header);

		//    if( path != "." )
		pop(xml);
	}
	
}

