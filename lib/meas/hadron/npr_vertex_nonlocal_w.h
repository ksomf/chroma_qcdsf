/*! \file npr_vertex_nonlocal_w.h
 *  \author Stanislav Kazmin (University Leipzig)
 *  \date 2015-11-24
 *  \brief NPR vertex calculations for the non local axial current.
 */

#ifndef __npr_vertex_nonlocal_w_h__
#define __npr_vertex_nonlocal_w_h__

#include "chromabase.h"

namespace Chroma
{
//! Used to Set Requested Link Patterns
/*! \ingroup hadron */
typedef void (*BBLinkPattern)(bool& DoThisPattern, bool& DoFurtherPatterns, multi1d< int >& LinkPattern);

//! NPR vertices
/*! \ingroup hadron */
void NprVertexNonlocal(const LatticePropagator& F,
					   const multi1d<LatticeColorMatrix>& U,
					   const unsigned short int MaxNLinks,
					   const BBLinkPattern LinkPattern,
					   QDPFileWriter& qio_file);

}  // end namespace Chroma

#endif // __npr_vertex_nonlocal_w_h__

