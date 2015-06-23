// -*- C++ -*-
// $Id: qprop_qcdsf_io.h, v 1.0 2011-08-31 15:26:29 bglaessle Exp $
/*! \file
 * \brief Modified structs for source / sink Metadata, for resmeared objects
 */

#ifndef __qprop_qcdsf_io_h__
#define __qprop_qcdsf_io_h__

#include "chromabase.h"
#include "io/qprop_io.h"

namespace Chroma {

	// HACK !!!!
	// "extern" somewhere in io/qprop_io.cc
	// but NOT in io/qprop_io.h !!!
	void readGaugeHeader(XMLReader& paramtop, const std::string& path, std::string& gauge_header);

	struct ResmearSourcePropQCDSF_t {
	//	PropSourceConst_t			source_header; // not necessary since a new one is constructed automatically
		PropSourceConst_t			orig_source_header;
		GroupXML_t					fake_smear_header;
		multi1d<GroupXML_t>			smearing_headers;

		std::string					gauge_header;
	};

	struct ResmearForwardPropQCDSF_t {
	//	PropSinkSmear_t		sink_header; // not necessary since a new one is constructed automatically
		PropSinkSmear_t		orig_sink_header;
		GroupXML_t			fake_smear_header;
		multi1d<GroupXML_t>	smearing_headers;

		ChromaProp_t		prop_header;
		PropSourceConst_t	source_header;
		std::string			gauge_header;
	};

	/*
	struct LinCombinationSourcePropQCDSF_t { 	// not used yet!!
		PropSourceConst_t	source_header;
		std::string			gauge_header;
	};

	struct LinCombinationSourcePropQCDSF_t { 	// not used yet!!
		ChromaProp_t		prop_header;
		PropSourceConst_t	source_header;
		PropSinkSmear_t		sink_header;
		std::string			gauge_header;
	};
	*/

	void readSmearingArray(XMLReader& xml, const std::string& path, multi1d<GroupXML_t>& arr );
	void read(XMLReader& xml, const std::string& path, ResmearSourcePropQCDSF_t& param);
	void read(XMLReader& xml, const std::string& path, ResmearForwardPropQCDSF_t& param);

	void writeSmearingArray(XMLWriter& xml, const std::string& path, const multi1d<GroupXML_t>& arr );
	void write(XMLWriter& xml, const std::string& path, const ResmearSourcePropQCDSF_t& param);
	void write(XMLWriter& xml, const std::string& path, const ResmearForwardPropQCDSF_t& param);
}

#endif
