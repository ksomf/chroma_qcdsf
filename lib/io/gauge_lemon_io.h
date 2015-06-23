// $Id: gauge_lemon_io.h,v 0.9 2011-05-15 04:58:55 bglaessle Exp $

/*! \file
 *  \brief Gauge reader/writers in ILDG format with Lemon
 */


#ifndef __gauge_lemon_io_h__
#define __gauge_lemon_io_h__

#include "chromabase.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

    struct ILDGFormat_t {
        string version;
        string field;
        int precision;
        int lx, ly, lz, lt;
    
        ILDGFormat_t() : precision(-2) {}
    };
    
    void read( XMLReader& paramtop, const std::string& path, ILDGFormat_t& ildgf );
    
    
    //! Read a gauge config in QIO format
    /*!
     * \ingroup io
     *
     * \param file_xml     xml reader holding config info ( Modify )
     * \param record_xml   xml reader holding config info ( Modify )
     * \param u            gauge configuration ( Modify )
     * \param file         path ( Read )
     * \param serpar       either QDPIO_SERIAL, QDPIO_PARALLEL ( Read )
     */    
    
    void readGaugeILDGLemon(XMLReader& file_xml, 
    		   XMLReader& record_xml,
    		   multi1d<LatticeColorMatrix>& u, 
    		   const std::string& file, 
    		   QDP_serialparallel_t serpar);

}  // end namespace Chroma

#endif // QDP_USE_LEMON

#endif // __gauge_lemon_io_h__

