// $Id: wuppertal_source_quark_smearing.cc, v 1.0 2011-08-31 16:19:59 bglaessle $
/*! \file
 *  \brief Wuppertal SOURCE smearing "operator", optimizes:
 		- if source(propagator) is spin-diagonal, only 1 of 16 spin components is smeared
 		- if source is on a single timeslice, only that timeslice is smeared
 */

#include "chromabase.h"

#include "meas/smear/quark_smearing_factory.h"
#include "meas/smear/wuppertal_source_quark_smearing.h"
#include "meas/smear/wuppertal_source_smear.h"

namespace QDP {

	Set TimeSliceSetQCDSF;

	struct TimeSliceFunc : public SetFunc {
		TimeSliceFunc(int dir): dir_decay(dir) {}

		int operator() (const multi1d<int>& coordinate) const {
			if ( (dir_decay<0) or (dir_decay>=Nd) ) {
				return 0 ;
			} else {
				return coordinate[dir_decay];
			}
		}

		int numSubsets() const {
			if ( (dir_decay<0) or (dir_decay>=Nd) ) {
				return 1 ;
			} else {
				return Layout::lattSize()[dir_decay] ;
			}
		}

		int dir_decay; // state
	};
	
	void initTimeSliceSet() {
		TimeSliceSetQCDSF.make( TimeSliceFunc(Nd-1) );
	}
}

namespace Chroma 
{
  // Read parameters
  void read (XMLReader& xml, const string& path, WuppertalSourceQuarkSmearingEnv::Params& param)
  {
    WuppertalSourceQuarkSmearingEnv::Params tmp (xml, path);
    param = tmp;
  }

  //! Parameters for running code
  void write (XMLWriter& xml, const string& path, const WuppertalSourceQuarkSmearingEnv::Params& param)
  {
    param.writeXML (xml, path);
  }


  //! Hooks to register the class
  namespace WuppertalSourceQuarkSmearingEnv
  {
    //! Callback function
    QuarkSmearing<LatticePropagator>* createProp (XMLReader& xml_in,
						  const std::string& path)
    {
      return new QuarkSmear<LatticePropagator> (Params (xml_in, path));
    }

    //! Callback function
    QuarkSmearing<LatticeStaggeredPropagator>* createStagProp (XMLReader& xml_in,
							       const std::string& path)
    {
      return new QuarkSmear<LatticeStaggeredPropagator> (Params (xml_in, path));
    }

    //! Callback function
    QuarkSmearing<LatticeFermion>* createFerm (XMLReader& xml_in,
					       const std::string& path)
    {
      return new QuarkSmear<LatticeFermion> (Params (xml_in, path));
    }
    
    //! Callback function
    QuarkSmearing<LatticeColorVector>* createColorVec (XMLReader& xml_in,
						       const std::string& path)
    {
      return new QuarkSmear<LatticeColorVector> (Params (xml_in, path));
    }
    
    //! Name to be used
    const std::string name = "WUPPERTAL_SOURCE_SMEAR";

    //! Local registration flag
    static bool registered = false;
    static bool registeredQDPGlobals = false;

    //! Register all the factories
    bool registerAll() 
    {
      bool success = true; 
      if (! registered)
      {
	success &= Chroma::ThePropSmearingFactory::Instance ().registerObject (name, createProp);
	success &= Chroma::TheStagPropSmearingFactory::Instance ().registerObject (name, createStagProp);
	success &= Chroma::TheFermSmearingFactory::Instance ().registerObject (name, createFerm);
	success &= Chroma::TheColorVecSmearingFactory::Instance ().registerObject (name, createColorVec);
	registered = true;
      }
      return success;
    }

    //! Register all the factories
    void registerQDPGlobals() 
    {
      if ( !registeredQDPGlobals ) {
		QDP::initTimeSliceSet();
		registeredQDPGlobals = true;
      }
    }


    //! Parameters for running code
    Params::Params (XMLReader& xml, const string& path)
    {
      XMLReader paramtop (xml, path);

      read (paramtop, "wvf_param", wvf_param);
      read (paramtop, "wvfIntPar", wvfIntPar);
      read (paramtop, "no_smear_dir", no_smear_dir);
      read (paramtop, "t_source", t_source );
      read (paramtop, "smearSingleSpin", smearSingleSpin );
      read (paramtop, "smearSingleTimeslice", smearSingleTimeslice );
    }


    //! Parameters for running code
    void Params::writeXML (XMLWriter& xml, const string& path) const
    {
      push (xml, path);
    
      write (xml, "wvf_kind", WuppertalSourceQuarkSmearingEnv::name);
      write (xml, "wvf_param", wvf_param);
      write (xml, "wvfIntPar", wvfIntPar);
      write (xml, "no_smear_dir", no_smear_dir);
      write (xml, "t_source", t_source );
      write (xml, "smearSingleTimeslice", smearSingleTimeslice );
      write (xml, "smearSingleSpin", smearSingleSpin );

      pop (xml);
    }


    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticePropagator>::operator() (LatticePropagator& quark,
					       const multi1d<LatticeColorMatrix>& u) const
    {
	  if( params.smearSingleSpin )
		  wuppertalSourceSmear_spinDiagonal (u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir, tslice );
	  else
	      wuppertalSourceSmear (u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir, tslice );
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeStaggeredPropagator>::operator() (LatticeStaggeredPropagator& quark,
						        const multi1d<LatticeColorMatrix>& u) const
    {
      wuppertalSourceSmear (u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir, tslice );
    }

    //! Smear the quark
    template<>
    void
    QuarkSmear<LatticeFermion>::operator() (LatticeFermion& quark,
					    const multi1d<LatticeColorMatrix>& u) const
    {
      wuppertalSourceSmear (u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir, tslice );
    }

    //! Smear the color-vector
    template<>
    void
    QuarkSmear<LatticeColorVector>::operator() (LatticeColorVector& quark,
					        const multi1d<LatticeColorMatrix>& u) const
    {
      wuppertalSourceSmear (u, quark, params.wvf_param, params.wvfIntPar, params.no_smear_dir, tslice );
    }

  }  // end namespace
}  // end namespace Chroma

