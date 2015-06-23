// $Id: propagator_lemon_io.cc,v 0.9 2011-05-15 04:58:55 bglaessle Exp $
/*! \file
 * \brief Routines associated with Chroma propagator IO
 */

#include "chromabase.h"
#include "io/propagator_lemon_io.h"
#include "io/lemon_util.h"
#include "util/ferm/transf.h"
#include "dml.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

template<typename T>
void readILDGBinaryData( QLemonReader& r, LatticePropagator& p, DML_Checksum* this_checksum ) {
    multi1d<int> subLattSize = Layout::subgridLattSize(); // "node" latt size
    multi1d<int> fact(Nd); 
    fact[0]=1;
    fact[1]=subLattSize[0];
    fact[2]=fact[1]*subLattSize[1];
    fact[3]=fact[2]*subLattSize[2];

    // dml stuff
    int latsize[4];
    latsize[0] = Layout::lattSize()[0];
    latsize[1] = Layout::lattSize()[1];
    latsize[2] = Layout::lattSize()[2];
    latsize[3] = Layout::lattSize()[3];
    
    DML_checksum_init( this_checksum );
    
    // size info
    uint64_t siteSize = Nd*Nd*Nc*Nc*2; 
    uint64_t precSize = sizeof(T);
        
    int totSites = Layout::sitesOnNode();
    int whoami = Layout::nodeNumber();

    QDPIO::cout << "LEMON: reading " << precSize
                << "-byte floatingpoint numbers" << endl;

    T* buff = new T[totSites*siteSize];
    
    int mapping[] = { 3, 2, 1, 0 };
    r.readParallelMapped( (void*) buff, siteSize*precSize, mapping );
    // r.readParallel( (void*) buff, siteSize*precSize );


    for(int i=0 ; i<totSites; i++) {
        multi1d<int> site = Layout::siteCoords(whoami,i);

        if( whoami == Layout::nodeNumber(site) ) {
            int coords[4];
            
            uint64_t offset = 0;
            for( int mu=0; mu<Nd; ++mu ) {
                int local_coord_mu = site[mu] % subLattSize[mu];
                offset += siteSize*fact[mu]*local_coord_mu;

                coords[mu] = site[mu];
            }
            T* ptr = buff + offset;

//          DML_SiteRank DML_lex_rank(const int coords[], int latdim, int latsize[]);
//          void DML_checksum_init(DML_Checksum *checksum)
//          void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank,
//           961                         char *buf, size_t size){
//           void DML_checksum_combine(DML_Checksum *checksum){
            {
                DML_SiteRank dml_siterank = DML_lex_rank( coords, Nd, latsize );
                DML_checksum_accum( this_checksum, dml_siterank, (char*) ptr, siteSize*precSize );
            }
            
            Propagator prop;
            ColorMatrix cm;
            for( int sA=0; sA<Nd; ++sA ) {
            for( int sB=0; sB<Nd; ++sB ) {

                for( int cA=0; cA<Nc; ++cA ) {
                for( int cB=0; cB<Nc; ++cB ) {
                    T re = fixEndianess( ptr[72*sA+18*sB+6*cA+2*cB] );
                    T im = fixEndianess( ptr[72*sA+18*sB+6*cA+2*cB+1] );
            
                    Complex val = cmplx(Real(re),Real(im));
                    pokeColor ( cm, val, cA, cB );
                } } // for cB, for cA
            
                pokeSpin( prop, cm, sA, sB );
            } } // for sB, sA

            //pokeSite(p,prop,site);
            p.elem(i) = prop.elem();
        } // if ( site @ whoami ) 
    } // for i
        
    delete [] buff;
    
    DML_checksum_combine( this_checksum );
}

void readLemonQIOLatProp(XMLReader& file_xml,
		   XMLReader& record_xml, 
		   LatticePropagator& p, 
		   const string& file ) 
{
    // well defined state of p
    p = zero;
    
    DML_Checksum file_checksum;
    DML_Checksum obj_checksum;
    SCIDACFileFormat_t file_format;
    SCIDACRecordFormat_t record_format;

    bool have_privatefile_hdr = false;
    bool have_file_hdr = false;
    bool have_privaterecord_hdr = false;
    bool have_record_hdr = false;
    bool have_data = false;
    bool have_checksum = false;

    QLemonReader r( file.c_str() );

    do {
        string recordname = r.recordName();
        QDPIO::cout << "record '" << recordname << "': len=" << r.recordSize() << endl;

        if ( recordname == "scidac-private-file-xml" ) {
            try {
                istringstream strm( readRecordToString(r) );
                //QDPIO::cout << "AAA\n" << strm.str() << "\nOOO" << endl;

                XMLReader xml( strm );
                read( xml, "/scidacFile", file_format );
                
                if ( file_format.volfmt != 0 ) {
                    QDPIO::cerr << "not SINGLE_FILE" << endl;
                    QDP_abort( 1 );
                }
                
                if ( file_format.spacetime != Nd ) {
                    QDPIO::cerr << "wrong number of dims " << file_format.spacetime << endl;
                    QDP_abort( 1 );
                }

                for ( int mu=0; mu<Nd; ++mu )
                    if ( file_format.dims[mu] != Layout::lattSize()[mu] ) {
                    QDPIO::cerr << "invalid size of dimension " << mu << ": " << file_format.dims[mu] << endl;
                    QDP_abort( 1 );
                }
            }
            catch ( const string& err ) {
                QDPIO::cerr << "caught exception in "
                            << __FILE__ << "@l" << __LINE__
                            << ": " << err << endl;
            
                //QDPIO::cerr << "proceeding (unsafe)" << endl;
                QDP_abort( 1 );
            } 
            have_privatefile_hdr=true;
        }
        else if ( recordname == "scidac-file-xml" ) {
            istringstream strm( readRecordToString(r) );
            file_xml.open( strm );
            have_file_hdr=true;
        }
        else if ( recordname == "scidac-private-record-xml" ) {
            try {
                istringstream strm( readRecordToString(r) );
                //QDPIO::cout << "AAA\n" << strm.str() << "\nOOO" << endl;
                
                XMLReader xml( strm );
                read( xml, "/scidacRecord", record_format );

                if ( record_format.datatype != "Lattice" ) {
                    QDPIO::cerr << "invalid datatype: "
                                << record_format.datatype << endl;
                    QDP_abort( 1 );
                }
                else if ( record_format.precision != "F" 
                       and record_format.precision != "D" ) {
                    QDPIO::cerr << "unknown precison: "
                                << record_format.precision << endl;
                    QDP_abort( 1 );
                }
                else if ( record_format.colors != Nc ) {
                    QDPIO::cerr << "wrong number of colors: "
                                << record_format.colors << endl;
                    QDP_abort( 1 );
                }
                else if ( record_format.spins != Ns ) {
                    QDPIO::cerr << "wrong number of spins: "
                                << record_format.spins << endl;
                    QDP_abort( 1 );
                }
                else if ( record_format.datacount != 1 ) {
                    QDPIO::cerr << "expected datacount 1 instead of "
                                << record_format.datacount << endl;
                    QDP_abort( 1 );
                }
            }
            catch ( const string& err ) {
                QDPIO::cerr << "caught exception in "
                            << __FILE__ << "@l" << __LINE__
                            << ": " << err << endl;
            
                //QDPIO::cerr << "proceeding (unsafe)" << endl;
                QDP_abort( 1 );
            } 
            have_privaterecord_hdr=true;
        }
        else if ( recordname == "scidac-record-xml" ) {
            istringstream strm( readRecordToString(r) );
            record_xml.open( strm );
            have_record_hdr=true;
        }
        else if ( recordname == "scidac-binary-data" ) {
            if ( have_privatefile_hdr and have_privaterecord_hdr ) {
                if( r.recordSize() != Layout::vol()*record_format.typesize ) {
                    QDPIO::cerr << "invalid recordsize " << r.recordSize() << endl;
                    QDP_abort( 1 );
                }

                if ( record_format.precision=="F" ) {
                    if ( record_format.typesize != Nc*Nc*Ns*Ns*8 ) {
                        QDPIO::cerr << "invalid typesize "
                                    << record_format.typesize << endl;
                        QDP_abort( 1 );
                    }

                    readILDGBinaryData<float>( r, p, &obj_checksum );
                }
                else if ( record_format.precision=="D" ) {
                    if ( record_format.typesize != Nc*Nc*Ns*Ns*16 ) {
                        QDPIO::cerr << "invalid typesize "
                                    << record_format.typesize << endl;
                        QDP_abort( 1 );
                    }

                    readILDGBinaryData<double>( r, p, &obj_checksum );
                }
            }
            else {
                QDPIO::cerr << __FILE__ << "@l" << __LINE__
                            << ": no previous file/record records"
                            << endl;
                QDP_abort(1);
            }
            have_data=true;
        }
        else if ( recordname == "scidac-checksum" ) {
            try {
                istringstream strm( readRecordToString(r) );
                //QDPIO::cout << "AAA\n" << strm.str() << "\nOOO" << endl;

                XMLReader xml( strm );
                read( xml, "/scidacChecksum", file_checksum );
            }
            catch ( const string& err ) {
                QDPIO::cerr << "caught exception in "
                            << __FILE__ << "@l" << __LINE__
                            << ": " << err << endl;
            
                //QDPIO::cerr << "proceeding (unsafe)" << endl;
                QDP_abort( 1 );
            } 
            have_checksum=true;
        }

        if( have_privatefile_hdr and
            have_file_hdr and
            have_privaterecord_hdr and
            have_record_hdr and
            have_data and
            have_checksum )
            break;
    
        r.nextRecord(); // will abort programm if not successful = no more records
    
    } while ( true );

    if ( file_checksum.suma != obj_checksum.suma
         or file_checksum.sumb != obj_checksum.sumb ) {

        QDPIO::cerr << __func__ << " - checksums dont match:\nFileChckSum "
                << file_checksum.suma << " | "
                << file_checksum.sumb << "\nObjChckSum  "
                << obj_checksum.suma << " | "
                << obj_checksum.sumb << endl;
        
        QDP_abort( 1 );
    }
}

void writeQIOBinaryData( QLemonWriter& w, const LatticePropagator& p, DML_Checksum* this_checksum ) {
    multi1d<int> subLattSize = Layout::subgridLattSize(); // "node" latt size
    multi1d<int> fact(Nd); 
    fact[0]=1;
    fact[1]=subLattSize[0];
    fact[2]=fact[1]*subLattSize[1];
    fact[3]=fact[2]*subLattSize[2];

    // dml stuff
    int latsize[4];
    latsize[0] = Layout::lattSize()[0];
    latsize[1] = Layout::lattSize()[1];
    latsize[2] = Layout::lattSize()[2];
    latsize[3] = Layout::lattSize()[3];
    
    DML_checksum_init( this_checksum );
    
    // size info
    uint64_t siteSize = Nd*Nd*Nc*Nc*2; 
    uint64_t precSize = sizeof(REAL);
        
    int totSites = Layout::sitesOnNode();
    int whoami = Layout::nodeNumber();

    QDPIO::cout << "LEMON: writing " << precSize
                << "-byte floatingpoint numbers" << endl;

    REAL* buff = new REAL[totSites*siteSize];
    
    for(int i=0 ; i<totSites; i++) {
        multi1d<int> site = Layout::siteCoords(whoami,i);

        if( whoami == Layout::nodeNumber(site) ) {
            int coords[4];
            
            uint64_t offset = 0;
            for( int mu=0; mu<Nd; ++mu ) {
                int local_coord_mu = site[mu] % subLattSize[mu];
                offset += siteSize*fact[mu]*local_coord_mu;

                coords[mu] = site[mu];
            }
            REAL* ptr = buff + offset;

            // Propagator prop = peekSite(p,site); // don't do this, avoid unnecessary mpi calls
            Propagator prop;
            prop.elem() = p.elem( i );

            for( int sA=0; sA<Nd; ++sA ) {
            for( int sB=0; sB<Nd; ++sB ) {
                ColorMatrix cm = peekSpin( prop, sA, sB );    

                for( int cA=0; cA<Nc; ++cA ) {
                for( int cB=0; cB<Nc; ++cB ) {
                    Complex val = peekColor( cm, cA, cB );
                    REAL raw_re = val.elem().elem().elem().real();
                    REAL raw_im = val.elem().elem().elem().imag();

                    ptr[72*sA+18*sB+6*cA+2*cB] = fixEndianess( raw_re );
                    ptr[72*sA+18*sB+6*cA+2*cB+1] = fixEndianess( raw_im );
                } } // for cB, for cA
            } } // for sB, sA

//          DML_SiteRank DML_lex_rank(const int coords[], int latdim, int latsize[]);
//          void DML_checksum_init(DML_Checksum *checksum)
//          void DML_checksum_accum(DML_Checksum *checksum, DML_SiteRank rank,
//           961                         char *buf, size_t size){
//           void DML_checksum_combine(DML_Checksum *checksum){
            DML_SiteRank dml_siterank = DML_lex_rank( coords, Nd, latsize );
            DML_checksum_accum( this_checksum, dml_siterank, (char*) ptr, siteSize*precSize );

        } // if ( site @ whoami ) 
    } // for i
      
    DML_checksum_combine( this_checksum );

    uint64_t msg_size = siteSize*precSize;
    msg_size *= Layout::vol();
    QDPIO::cout << "writing record 'scidac-binary-data': len=" << msg_size << endl;
    w.setRecordHeader( "scidac-binary-data", msg_size, 0, 0 );

    int mapping[] = { 3, 2, 1, 0 };
    w.writeParallelMapped( (void*) buff, siteSize*precSize, mapping );
    w.endRecord();
    // w.writeParallel( (void*) buff, siteSize*precSize );

    delete [] buff;
}

void writeLemonQIOLatProp(XMLBufferWriter& file_xml,
		   XMLBufferWriter& record_xml, 
		   const LatticePropagator& p, 
		   const string& file ) 
{
    DML_Checksum the_checksum;

    QLemonWriter w( file.c_str() );

    uint64_t msg_size;

    {
        SCIDACFileFormat_t file_format;
        file_format.version="1.1";
        file_format.spacetime=Nd;
        file_format.dims=Layout::lattSize();
        file_format.volfmt=0; // SINGLE_FILE

        XMLBufferWriter xmlbuff;
        write( xmlbuff, "scidacFile", file_format );
        msg_size=xmlbuff.str().size();
        MPI_Bcast( (char*)&msg_size, 8, MPI_CHAR, 0, QLEMON_COMM_CART);

        QDPIO::cout << "writing record 'scidac-private-file-xml': len=" << msg_size << endl;
        w.setRecordHeader( "scidac-private-file-xml", msg_size, 1, 0 );
        w.writeSerial( (void*)xmlbuff.str().c_str(), msg_size );
        w.endRecord();
    }

    {
        msg_size=file_xml.str().size();
        MPI_Bcast( (char*)&msg_size, 8, MPI_CHAR, 0, QLEMON_COMM_CART);

        QDPIO::cout << "writing record 'scidac-file-xml': len=" << msg_size << endl;
        w.setRecordHeader( "scidac-file-xml", msg_size, 0, 1 );
        w.writeSerial( (void*)file_xml.str().c_str(), msg_size );
        w.endRecord();
    }
    
    {
        SCIDACRecordFormat_t record_format;
        record_format.version="1.1";
        
        time_t cu_time;
        time(&cu_time);
        record_format.date=asctime(gmtime(&cu_time));
        size_t n = record_format.date.size();
        if ( record_format.date[n-1]=='\n' )
            record_format.date[n-1]=' ';
        record_format.date+="UTC";

        record_format.recordtype=0; // full lattice object
        record_format.datatype="Lattice";
        record_format.precision= sizeof(REAL)==8 ? "D" : "F";
        record_format.colors=Nc;
        record_format.spins=Ns;
        record_format.typesize=Nc*Nc*Ns*Ns*2*sizeof(REAL);
        record_format.datacount=1;

        XMLBufferWriter xmlbuff;
        write( xmlbuff, "scidacRecord", record_format );
        msg_size=xmlbuff.str().size();
        MPI_Bcast( (char*)&msg_size, 8, MPI_CHAR, 0, QLEMON_COMM_CART);

        QDPIO::cout << "writing record 'scidac-private-record-xml': len=" << msg_size << endl;
        w.setRecordHeader( "scidac-privte-record-xml", msg_size, 1, 0 );
        w.writeSerial( (void*)xmlbuff.str().c_str(), msg_size );
        w.endRecord();
    }

    {
        msg_size=record_xml.str().size();
        MPI_Bcast( (char*)&msg_size, 8, MPI_CHAR, 0, QLEMON_COMM_CART);

        QDPIO::cout << "writing record 'scidac-record-xml': len=" << msg_size << endl;
        w.setRecordHeader( "scidac-record-xml", msg_size, 0, 0 );
        w.writeSerial( (void*)record_xml.str().c_str(), msg_size );
        w.endRecord();
    }

    writeQIOBinaryData( w, p, &the_checksum );

    {
        XMLBufferWriter xmlbuff;
        write( xmlbuff, "scidacChecksum", the_checksum );
        msg_size=xmlbuff.str().size();
        MPI_Bcast( (char*)&msg_size, 8, MPI_CHAR, 0, QLEMON_COMM_CART);

        QDPIO::cout << "writing record 'scidac-checksum': len=" << msg_size << endl;
        w.setRecordHeader( "scidac-checksum", msg_size, 0, 1 );
        w.writeSerial( (void*)xmlbuff.str().c_str(), msg_size );
        w.endRecord();
    }
}

}  // end namespace Chroma

#endif // QDP_USE_LEMON

