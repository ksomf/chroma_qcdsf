// $Id: propagator_lemon_io.h,v 0.9 2011-05-15 04:58:55 bglaessle Exp $

/*! \file
 *  \brief Propagator reader/writers in QIO/ILDG format with Lemon
 */


#ifndef __propagator_lemon_io_h__
#define __propagator_lemon_io_h__

#include "chromabase.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

    //! Read a propagator in QIO format
    /*!
     * \ingroup io
     *
     * \param file_xml     xml reader holding config info ( Modify )
     * \param record_xml   xml reader holding config info ( Modify )
     * \param p            propagator ( Modify )
     * \param file         path ( Read )
     */    
    
    void readLemonQIOLatProp(XMLReader& file_xml, 
    		   XMLReader& record_xml,
    		   LatticePropagator& p, 
    		   const std::string& file );

    void writeLemonQIOLatProp(XMLBufferWriter& file_xml, 
    		   XMLBufferWriter& record_xml,
    		   const LatticePropagator& p, 
    		   const std::string& file );

}  // end namespace Chroma

#endif // QDP_USE_LEMON

#endif // __gauge_lemon_io_h__

