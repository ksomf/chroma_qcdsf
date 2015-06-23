// $Id: inline_resmear_sink_qcdsf_w.cc,v 1.0 2011-08-31 12:54:23 bglaessle Exp $
/*! \file
 * \brief Inline construction of resmear_sink_qcdsf
 *
 * Sink smear propagators
 */

#include "handle.h"
#include "meas/inline/hadron/inline_resmear_sink_qcdsf_w.h"
#include "meas/inline/abs_inline_measurement_factory.h"
#include "meas/sinks/sink_smearing_factory.h"
#include "meas/sinks/sink_smearing_aggregate.h"
#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/quark_smearing_aggregate.h"
#include "meas/glue/mesplq.h"
#include "util/ft/sftmom.h"
#include "util/info/unique_id.h"
#include "io/qprop_qcdsf_io.h"

#include "meas/inline/io/named_objmap.h"

//#define DEBUG_INLINE_RESMEAR_SINK_QCSF_W

namespace Chroma {

    //! Propagator input
    void read( XMLReader& xml, const string& path, InlineResmearSinkEnvQCDSF::Params::NamedObject_t& input ) {
        XMLReader inputtop( xml, path );

        read( inputtop, "smeared_gauge_id", input.gauge_id );
        read( inputtop, "input_sink_id", input.input_prop_id );
        read( inputtop, "output_sink_id", input.output_prop_id );
    }

    //! Propagator output
    void write( XMLWriter& xml, const string& path, const InlineResmearSinkEnvQCDSF::Params::NamedObject_t& input ) {
        push( xml, path );

        write( xml, "smeared_gauge_id", input.gauge_id );
        write( xml, "input_sink_id", input.input_prop_id );
        write( xml, "output_sink_id", input.output_prop_id );

        pop( xml );
    }

    //! Propagator input
    void read( XMLReader& xml, const string& path, InlineResmearSinkEnvQCDSF::Params::ResmearParam_t& input ) {
        XMLReader inputtop( xml, path );

        read( inputtop, "version", input.version );

        switch ( input.version ) {
                /**************************************************************************/
            case 1:
                read( inputtop, "j_decay", input.j_decay );
                input.quark_smearing = readXMLGroup( inputtop, "SmearingParam", "wvf_kind" );
                break;

                /**************************************************************************/
            default:
                QDPIO::cerr << "Resmear parameter version " << input.version << " unsupported." << endl;
                QDP_abort( 1 );
        }
    }

    //! Propagator output
    void write( XMLWriter& xml, const string& path, const InlineResmearSinkEnvQCDSF::Params::ResmearParam_t& input ) {
        push( xml, path );

        int version = 1;
        write( xml, "version", version );
        write( xml, "j_decay", input.j_decay );
        xml << input.quark_smearing.xml;

        pop( xml );
    }


    namespace InlineResmearSinkEnvQCDSF {
        namespace {
            AbsInlineMeasurement* createMeasurement( XMLReader& xml_in,
                    const std::string& path ) {
                return new InlineMeas( Params( xml_in, path ) );
            }

            //! Local registration flag
            bool registered = false;
        }

        const std::string name = "RESMEAR_SINK-QCDSF";

        //! Register all the factories
        bool registerAll() {
            bool success = true;

            if ( ! registered ) {
                success &= QuarkSinkSmearingEnv::registerAll();
                success &= TheInlineMeasurementFactory::Instance().registerObject( name, createMeasurement );
                registered = true;
            }

            return success;
        }


        // Param stuff
        Params::Params() {
            frequency = 0;
        }

        Params::Params( XMLReader& xml_in, const std::string& path ) {
            try {
                XMLReader paramtop( xml_in, path );

                if ( paramtop.count( "Frequency" ) == 1 )
                    read( paramtop, "Frequency", frequency );
                else
                    frequency = 1;

                // Parameters for source construction
                read( paramtop, "Param", param );

                // Read in the output propagator/source configuration info
                read( paramtop, "NamedObject", named_obj );
            }
            catch ( const std::string& e ) {
                QDPIO::cerr << __func__ << ": Caught Exception reading XML: " << e << endl;
                QDP_abort( 1 );
            }
        }


        void
        Params::writeXML( XMLWriter& xml_out, const std::string& path ) {
            push( xml_out, path );

            // Parameters for source construction
            write( xml_out, "Param", param );

            // Write out the output propagator/source configuration info
            write( xml_out, "NamedObject", named_obj );

            pop( xml_out );
        }


        void
        InlineMeas::operator()( unsigned long update_no,
                                XMLWriter& xml_out ) {
            START_CODE();

            QDPIO::cout << name << ": resmearing for (propagator) sinks" << endl;

            StopWatch snoop;
            snoop.reset();
            snoop.start();

            // Test and grab a reference to the gauge field
            XMLBufferWriter gauge_xml;

            try {
                TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >( params.named_obj.gauge_id );
                TheNamedObjMap::Instance().get( params.named_obj.gauge_id ).getRecordXML( gauge_xml );
            }
            catch ( std::bad_cast ) {
                QDPIO::cerr << name << ": caught dynamic cast error"
                            << endl;
                QDP_abort( 1 );
            }
            catch ( const string& e ) {
                QDPIO::cerr << name << ": map call failed: " << e
                            << endl;
                QDP_abort( 1 );
            }

            const multi1d<LatticeColorMatrix>& u_smr =
                TheNamedObjMap::Instance().getData< multi1d<LatticeColorMatrix> >( params.named_obj.gauge_id );

            push( xml_out, "resmear_sink_qcdsf" );

            write( xml_out, "update_no", update_no );

            // Write out the input
            params.writeXML( xml_out, "Input" );

            // Write out the config header
            write( xml_out, "Config_info", gauge_xml );

            // Calculate some gauge invariant observables just for info.
            MesPlq( xml_out, "Observables", u_smr );

            //
            // Read the quark propagator and extract headers
            //
            XMLReader prop_file_xml, prop_record_xml;

            LatticePropagator quark_propagator;

            int j_decay;   // need this for diagnostics

            QDPIO::cout << "Attempt to read forward propagator" << endl;

            try {
                // Grab a copy of the propagator. Will modify it later.
                quark_propagator =
                    TheNamedObjMap::Instance().getData<LatticePropagator>( params.named_obj.input_prop_id );

                // Snarf the prop info. This is will throw if the input_prop_id is not there
                TheNamedObjMap::Instance().get( params.named_obj.input_prop_id ).getFileXML( prop_file_xml );
                TheNamedObjMap::Instance().get( params.named_obj.input_prop_id ).getRecordXML( prop_record_xml );

                // Snarf out the first j_decay
                read( prop_record_xml, "/descendant::j_decay[1]", j_decay );

                // Write out the (input) sink header
                write( xml_out, "Sink_file_info", prop_file_xml );
                write( xml_out, "Sink_record_info", prop_record_xml );
            }
            catch ( std::bad_cast ) {
                QDPIO::cerr << name << ": caught dynamic cast error"
                            << endl;
                QDP_abort( 1 );
            }
            catch ( const string& e ) {
                QDPIO::cerr << name << ": error extracting prop_header: " << e << endl;
                QDP_abort( 1 );
            }

            // Sanity check - write out the norm2 of the forward prop in the j_decay direction
            // Use this for any possible verification
            {
                // Initialize the slow Fourier transform phases
                SftMom phases( 0, true, j_decay );

                multi1d<Double> forward_prop_corr = sumMulti( localNorm2( quark_propagator ),
                                                    phases.getSet() );

                push( xml_out, "Forward_prop_correlator" );
                write( xml_out, "forward_prop_corr", forward_prop_corr );
                pop( xml_out );
            }

            //
            // resmear the SINK
            //
            try {
                //  QDPIO::cout << "Smearing_xml = " << params.param.quark_smearing.xml << endl;
                std::istringstream  xml_s( params.param.quark_smearing.xml );
                XMLReader smeartop( xml_s );


                Handle< QuarkSmearing<LatticePropagator> >
                quarkSmearing( ThePropSmearingFactory::Instance().createObject( params.param.quark_smearing.id,
                               smeartop,
                               params.param.quark_smearing.path ) );

                ( *quarkSmearing )( quark_propagator, u_smr );

                /*          std::istringstream  xml_s(params.param.source.xml);
                            XMLReader  sourcetop(xml_s);
                            QDPIO::cout << "Smearing = " << params.param.source.id << endl;

                            Handle< QuarkSmearing<LatticePropagator> >
                                sourceSmearing(ThePropSourceSmearingFactory::Instance().createObject(params.param.source.id,
                                                                                                    sourcetop,
                                                                                                    params.param.source.path,
                                                                                                    u)); */
            }
            catch ( const std::string& e )    {
                QDPIO::cerr << name << ": Caught Exception Resmearing source: " << e << endl;
                QDP_abort( 1 );
            }

            // Sanity check - write out the propagator (pion) correlator in the Nd-1 direction
            {
                // Initialize the slow Fourier transform phases
                SftMom phases( 0, true, j_decay );

                multi1d<Double> prop_corr = sumMulti( localNorm2( quark_propagator ),
                                                      phases.getSet() );

                push( xml_out, "SinkSmearedProp_correlator" );
                write( xml_out, "sink_smeared_prop_corr", prop_corr );
                pop( xml_out );
            }

            // Save the propagator
            try {
                XMLBufferWriter file_xml;
                push( file_xml, "resmear_sink_qcdsf" );
                write( file_xml, "id", uniqueId() ); // NOTE: new ID form
                pop( file_xml );

                XMLBufferWriter record_xml;

                // Construct an appropriate record_xml based on the input type
                if ( prop_record_xml.count( "/SinkSmear" ) != 0 ) {
                    ResmearForwardPropQCDSF_t  orig_header;
                    read( prop_record_xml, "/SinkSmear", orig_header );

                    int this_smear = orig_header.smearing_headers.size();
                    ResmearForwardPropQCDSF_t new_header = orig_header;
                    new_header.smearing_headers.resize( this_smear + 1 );

                    for ( size_t j = 0; j < orig_header.smearing_headers.size(); ++j )
                        new_header.smearing_headers[j] = orig_header.smearing_headers[j];

                    new_header.smearing_headers[this_smear] = params.param.quark_smearing;

                    write( record_xml, "SinkSmear", new_header );
                    //      write(xml_out, "SinkSmear", new_header);
                }
                else {
                    throw std::string( "No appropriate header found" );
                }

                // Write the smeared prop xml info
                TheNamedObjMap::Instance().create<LatticePropagator>( params.named_obj.output_prop_id );
                TheNamedObjMap::Instance().getData<LatticePropagator>( params.named_obj.output_prop_id )
                = quark_propagator;
                TheNamedObjMap::Instance().get( params.named_obj.output_prop_id ).setFileXML( file_xml );
                TheNamedObjMap::Instance().get( params.named_obj.output_prop_id ).setRecordXML( record_xml );

                QDPIO::cout << "Sink successfully updated" << endl;
            }
            catch ( std::bad_cast ) {
                QDPIO::cerr << name << ": dynamic cast error"
                            << endl;
                QDP_abort( 1 );
            }
            catch ( const string& e ) {
                QDPIO::cerr << name << ": error message: " << e << endl;
                QDP_abort( 1 );
            }

            pop( xml_out ); // resmear_sink_qcdsf


            snoop.stop();
            QDPIO::cout << name << ": total time = "
                        << snoop.getTimeInSeconds()
                        << " secs" << endl;

            QDPIO::cout << name << ": ran successfully" << endl;

            END_CODE();
        }

    }

}

#ifdef DEBUG_INLINE_RESMEAR_SINK_QCSF_W
#undef DEBUG_INLINE_RESMEAR_SINK_QCSF_W
#endif

