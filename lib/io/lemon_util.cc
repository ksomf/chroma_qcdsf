// $Id: lemon_util.cc,v 0.9 2011-05-15 04:58:55 bglaessle Exp $
/*! \file
 * \brief Routines associated with Lemon IO
 */

#include "chromabase.h"
#include "io/lemon_util.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

void read( XMLReader& paramtop, const std::string& path, SCIDACFileFormat_t& scidacf ) {
    XMLReader xml( paramtop, path );
    read( xml, "version", scidacf.version );
    read( xml, "spacetime", scidacf.spacetime );
    read( xml, "dims", scidacf.dims );
    read( xml, "volfmt", scidacf.volfmt );
}

void read( XMLReader& paramtop, const std::string& path, SCIDACRecordFormat_t& scidacf ) {
    XMLReader xml( paramtop, path );
    read( xml, "version", scidacf.version );
    read( xml, "date", scidacf.date );
    read( xml, "recordtype", scidacf.recordtype );
    read( xml, "datatype", scidacf.datatype );
    read( xml, "precision", scidacf.precision );
    read( xml, "colors", scidacf.colors );
    read( xml, "spins", scidacf.spins );
    read( xml, "typesize", scidacf.typesize );
    read( xml, "datacount", scidacf.datacount );
}

void read( XMLReader& paramtop, const std::string& path, DML_Checksum& checksum ) {
    XMLReader xml( paramtop, path );

    string hex_a, hex_b;
    read( xml, "suma", hex_a );
    read( xml, "sumb", hex_b );

    stringstream strm( hex_a+" "+hex_b );
    strm >> std::hex >> checksum.suma;
    strm >> std::hex >> checksum.sumb;
    //QDPIO::cout << "suma=" << hex_a << " -> " << checksum.suma
    //            << "\nsumb=" << hex_b << " -> " << checksum.sumb << endl;
}

void write( XMLWriter& xml, const std::string& path, const SCIDACFileFormat_t& scidacf ) {
    push( xml, path );
    write( xml, "version", scidacf.version );
    write( xml, "spacetime", scidacf.spacetime );
    write( xml, "dims", scidacf.dims );
    write( xml, "volfmt", scidacf.volfmt );
    pop( xml );
}

void write( XMLWriter& xml, const std::string& path, const SCIDACRecordFormat_t& scidacf ) {
    push( xml, path );
    write( xml, "version", scidacf.version );
    write( xml, "date", scidacf.date );
    write( xml, "recordtype", scidacf.recordtype );
    write( xml, "datatype", scidacf.datatype );
    write( xml, "precision", scidacf.precision );
    write( xml, "colors", scidacf.colors );
    write( xml, "spins", scidacf.spins );
    write( xml, "typesize", scidacf.typesize );
    write( xml, "datacount", scidacf.datacount );
    pop( xml );
}

void write( XMLWriter& xml, const std::string& path, const DML_Checksum& checksum ) {
    push( xml, path );

    write( xml, "version", "1.0" );

    stringstream strm;
    strm << std::hex << checksum.suma;
    write( xml, "suma", strm.str() );
    //QDPIO::cout << "suma=" << strm.str() << " -> " << checksum.suma << endl;

    strm.str("");
    strm << std::hex << checksum.sumb;
    write( xml, "sumb", strm.str() );
    //QDPIO::cout << "suma=" << strm.str() << " -> " << checksum.sumb << endl;

    pop( xml );
}

std::string readRecordToString( QLemonReader& r ) {
    uint64_t bsize = r.recordSize();
    char* cbuff = new char[bsize+1];
    cbuff[bsize] = 0;
    r.readSerial( (void*)cbuff, bsize );
    string sbuff(cbuff);
    //string sbuff(cbuff,bsize);
    delete [] cbuff;
    return sbuff; 
}

}  // end namespace Chroma

#endif // QDP_USE_LEMON

