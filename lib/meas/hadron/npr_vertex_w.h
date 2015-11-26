// -*- C++ -*-
//  $Id: npr_vertex_w.h,v 1.5 2008-12-21 21:22:37 edwards Exp $
/*! \file
 *  \brief NPR vertex calculations
 */

#ifndef __npr_vertex_w_h__
#define __npr_vertex_w_h__

#include "chromabase.h"

namespace Chroma 
{
	//! Used to Set Requested Link Patterns
	/*! \ingroup hadron */
	typedef void (*BBLinkPattern)(bool &                          DoThisPattern,
				bool &                          DoFurtherPatterns,
				multi1d< int > & LinkPattern);

	//! NPR vertices
	/*! \ingroup hadron */
	void NprVertex(const LatticePropagator &             F,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const BBLinkPattern                   LinkPattern,
		 QDPFileWriter& qio_file);
	
	  QDP::StandardOutputStream& operator<<(QDP::StandardOutputStream& s, const multi1d<int>& d);
	  
	  void AddLinks(const LatticePropagator&  B,
		const LatticePropagator&  F,
		const multi1d< LatticeColorMatrix > & U,
		multi1d< int >&    LinkDirs,
		const int          MaxNLinks,
		BBLinkPattern      LinkPattern,
		const int          PreviousDir,
		const int          PreviousMu,
		QDPFileWriter&     qio_file,
		int&               GBB_NLinkPatterns);
		
		  void NprVertex(const LatticePropagator &             F,
		 const multi1d< LatticeColorMatrix > & U,
		 const unsigned short int              MaxNLinks,
		 const BBLinkPattern                   LinkPattern,
		 QDPFileWriter& qio_file);

}  // end namespace Chroma

#endif

//###################################################################################//
//###################################################################################//
