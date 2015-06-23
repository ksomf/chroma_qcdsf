// $Id: lemon_util.h,v 0.9 2011-05-15 04:58:55 bglaessle Exp $

/*! \file
 *  \brief utilities for parallel IO of QIO e.g. Propagators with Lemon
 */


#ifndef __lemon_util_h__
#define __lemon_util_h__

#include "chromabase.h"
#include "dml.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

    std::string readRecordToString( QLemonReader& r );

    struct SCIDACFileFormat_t {
        string version;
        int spacetime;
        multi1d<int> dims;
        int volfmt;
    };

    struct SCIDACRecordFormat_t {
        string version;
        string date;
        int recordtype;
        string datatype;
        string precision;
        int colors;
        int spins;
        int typesize;
        int datacount;
    };

    void read( XMLReader& paramtop, const std::string& path, SCIDACFileFormat_t& scidacf );
    void read( XMLReader& paramtop, const std::string& path, SCIDACRecordFormat_t& scidacf );
    void read( XMLReader& paramtop, const std::string& path, DML_Checksum& checksum );

    void write( XMLWriter& xml, const std::string& path, const SCIDACFileFormat_t& scidacf );
    void write( XMLWriter& xml, const std::string& path, const SCIDACRecordFormat_t& scidacf );
    void write( XMLWriter& xml, const std::string& path, const DML_Checksum& checksum );

    static int machineIsLittleEndian() {
        union {
            int  l;
            char c[sizeof(int)];
        } u;
        u.l = 1;

        return (u.c[sizeof(int)-1] == 1 ? 0 : 1);
    }

    template<class T>
    T swapEndianess( const T inType  ) {
        T retVal;
        char *typeToConvert = ( char* ) & inType;
        char *returnType = ( char* ) & retVal;
        size_t N = sizeof(T);

        // swap the bytes into a temporary buffer
        for( size_t j=0; j<N; ++j )
            returnType[j] = typeToConvert[N-1-j];

        return retVal;
    }

    template<class T>
    inline T fixEndianess( const T inType  ) {
        if( machineIsLittleEndian()==1 )
            return swapEndianess( inType );

        return inType;
    }

}  // end namespace Chroma

#endif // QDP_USE_LEMON

#endif // __gauge_lemon_io_h__

