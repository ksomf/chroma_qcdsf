// $Id: gauge_lemon_io.cc,v 0.9 2011-05-15 04:58:55 bglaessle Exp $
/*! \file
 * \brief Routines associated with Chroma propagator gauge IO
 */

#include "chromabase.h"
#include "io/gauge_lemon_io.h"
#include "io/lemon_util.h"
#include "util/gauge/reunit.h"

#ifdef QDP_USE_LEMON

namespace Chroma {

void read( XMLReader& paramtop, const std::string& path, ILDGFormat_t& ildgf ) {
    XMLReader xml( paramtop, path );
    read( xml, "version", ildgf.version );
    read( xml, "field", ildgf.field );
    read( xml, "precision", ildgf.precision );
    read( xml, "lx", ildgf.lx );
    read( xml, "ly", ildgf.ly );
    read( xml, "lz", ildgf.lz );
    read( xml, "lt", ildgf.lt );
/*    read( paramtop, "/ildgFormat/version", ildgf.version );
    read( paramtop, "/ildgFormat/field", ildgf.field );
    read( paramtop, "/ildgFormat/precision", ildgf.precision );
    read( paramtop, "/ildgFormat/lx", ildgf.lx );
    read( paramtop, "/ildgFormat/ly", ildgf.ly );
    read( paramtop, "/ildgFormat/lz", ildgf.lz );
    read( paramtop, "/ildgFormat/lt", ildgf.lt ); 
    read( paramtop, "/descendant::version[1]", ildgf.version );
    read( paramtop, "/descendant::field[1]", ildgf.field );
    read( paramtop, "/descendant::precision[1]", ildgf.precision );
    read( paramtop, "/descendant::lx[1]", ildgf.lx );
    read( paramtop, "/descendant::ly[1]", ildgf.ly );
    read( paramtop, "/descendant::lz[1]", ildgf.lz );
    read( paramtop, "/descendant::lt[1]", ildgf.lt ); */
}

template<typename T>
void readILDGBinaryData( QLemonReader& r, multi1d<LatticeColorMatrix>& u ) {
    multi1d<int> subLattSize = Layout::subgridLattSize(); // "node" latt size
    multi1d<int> fact(Nd); 
    fact[0]=1;
    fact[1]=subLattSize[0];
    fact[2]=fact[1]*subLattSize[1];
    fact[3]=fact[2]*subLattSize[2];

    // Nd*(Nc^2)*cmplx = 72 for SU(3) in 3+1 Dimensions
    uint64_t siteSize = Nd*Nc*Nc*2; 
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
            uint64_t offset = 0;
            for( int mu=0; mu<Nd; ++mu ) {
                int local_coord_mu = site[mu] % subLattSize[mu];
                offset += siteSize*fact[mu]*local_coord_mu;
            }
            T* ptr = buff + offset;
            
            for( int mu=0; mu<Nd; ++mu ) {
                ColorMatrix cm;
                for( int colA = 0; colA<Nc; ++colA )
                for( int colB = 0; colB<Nc; ++colB ) {
                    T re = fixEndianess( ptr[18*mu+6*colA+2*colB] );
                    T im = fixEndianess( ptr[18*mu+6*colA+2*colB+1] );
            
                    Complex val = cmplx(Real(re),Real(im));
                    pokeColor ( cm, val, colA, colB );
                } // for colB, for colA
            
                //pokeSite(u[mu],cm,site);
                u[mu].elem(i) = cm.elem();
            } // for mu
        } // if ( site @ whoami ) 
    } // for i
        
    delete [] buff;
}

void readGaugeILDGLemon(XMLReader& file_xml,
		   XMLReader& record_xml, 
		   multi1d<LatticeColorMatrix>& u, 
		   const string& file, 
		   QDP_serialparallel_t serpar)
{
    // well defined state of u
    u.resize( Nd );
    u = zero;
    
    ILDGFormat_t ildg_format;
    string ildg_lfn;

    bool have_hdr = false;
    bool have_data = false;
    bool have_lfn = false;

    QLemonReader r( file.c_str() );

    do {
        string recordname = r.recordName();
        QDPIO::cout << "record '" << recordname << "': len=" << r.recordSize() << endl;

        if ( recordname == "ildg-format" ) {
            try {
                istringstream sstr( readRecordToString( r ) );
                //QDPIO::cout << "AAA\n" << sstr.str() << "\nOOO" << endl;

                XMLReader xml( sstr );
                if ( xml.count( "/ildgFormat" ) != 0 ) {
                    read( xml, "/ildgFormat", ildg_format );
                
                    if( ildg_format.precision != 32 and
                        ildg_format.precision != 64 ) {
                        QDPIO::cerr << "invalid ILDG floatingpoint precision: "
                                    << ildg_format.precision << endl;
                        QDP_abort( 1 );
                    }
            
                    if ( ildg_format.lx != Layout::lattSize()[0] 
                        or ildg_format.ly != Layout::lattSize()[1]  
                        or ildg_format.lz != Layout::lattSize()[2]  
                        or ildg_format.lt != Layout::lattSize()[3] ) {
                        QDPIO::cerr << "size of dimensions don't match" << endl;
                        QDP_abort( 1 );
                    }
                }
                else {
                    // on some ensembles there are <ildgFormat ...> tags
                    // which the XMLReader cannot handle
                    QDPIO::cerr << "couldn't find valid '/ildgFormat' in 'ildg-format' record\n"
                                << "proceeding (unsafe)" << endl;
                
                    // the only information really used from the ILDGFormat
                    // is the precision, which needs to be derived from
                    // the binary record-size if reading the meta-data fails.
                    ildg_format.precision = -1;
                }
            }
            catch ( const string& err ) {
                // exception was probably thrown by the XMLReader,
                QDPIO::cerr << "caught exception in "
                            << __FILE__ << "@l" << __LINE__
                            << ": " << err << endl;
                QDP_abort( 1 );
            } 
            have_hdr=true;
        }
        else if ( recordname == "ildg-binary-data" ) {
            if ( have_hdr ) {
                if( ildg_format.precision == -1 ) {
                    uint64_t rsize=r.recordSize();
                    uint64_t latVol=Layout::vol();

                    int siteSize=Nd*Nc*Nc*2;

                    if ( rsize%( latVol*siteSize ) == 0 ) {
                        ildg_format.precision = 8*rsize/latVol/siteSize;
                    }
                    else {
                        QDPIO::cerr << "unable to determine precision" << endl;
                        QDP_abort( 1 ) ;
                    }
                }

                switch (ildg_format.precision) {
                    case 32:
                        readILDGBinaryData<float>( r, u );
                        break;
                        
                    case 64:
                        readILDGBinaryData<double>( r, u );
                        break;
                        
                    default:
                        QDPIO::cerr << "unknown floating point precision "
                                    << ildg_format.precision << endl;
                        QDP_abort( 1 );
                }

            //    if ( ildg_format.precision != BASE_PRECISION ) {
            //        QDPIO::cout << "file precision != build precision: reunitarizing" << endl;
            //        for( int mu=0; mu<Nd; ++mu )
            //            reunit( u[mu] );
            //    }
            }
            else {
                QDPIO::cerr << __FILE__ << "@l" << __LINE__
                            << ": no previous 'ildg-format' record"
                            << endl;
                QDP_abort(1);
            }
            have_data=true;
        }
        else if ( recordname == "ildg-data-lfn" ) {
            ildg_lfn = readRecordToString( r );
            have_lfn=true;
        }

        if( have_hdr and have_data and have_lfn )
            break;
    
        r.nextRecord(); // will abort programm if not successful = no more records
    
    } while ( true );
  
    {                
        XMLBufferWriter f_xml;
        push(f_xml, "gauge");
        write(f_xml, "id", int(0) );
        pop(f_xml);
        file_xml.open(f_xml);
    }
     
    {
        XMLBufferWriter r_xml;
        push(r_xml, "ILDG");
        write(r_xml,"LFN",ildg_lfn );
        pop(r_xml);
        record_xml.open(r_xml);
    }
}


}  // end namespace Chroma

#endif // QDP_USE_LEMON

