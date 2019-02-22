//  $Id: mesons2_w.cc,v 1.1 2006-07-10 19:53:37 edwards Exp $
//  $Log: mesons2_w.cc,v $
//  Revision 1.1  2006-07-10 19:53:37  edwards
//  A complex version.
//
//  Revision 3.0  2006/04/03 04:59:00  edwards
//  Major overhaul of fermion and gauge action interface. Basically,
//  all fermacts and gaugeacts now carry out  <T,P,Q>  template parameters. These are
//  the fermion type, the "P" - conjugate momenta, and "Q" - generalized coordinates
//  in the sense of Hamilton's equations. The fermbc's have been rationalized to never
//  be over multi1d<T>. The "createState" within the FermionAction is now fixed meaning
//  the "u" fields are now from the coordinate type. There are now "ConnectState" that
//  derive into FermState<T,P,Q> and GaugeState<P,Q>.
//
//  Revision 2.1  2005/11/20 18:29:00  edwards
//  Added parenthesis
//
//  Revision 2.0  2005/09/25 21:04:35  edwards
//  Moved to version 2.0
//
//  Revision 1.21  2005/02/28 03:35:25  edwards
//  Fixed doxygen comments.
//
//  Revision 1.20  2005/01/14 18:42:36  edwards
//  Converted all lib files to be in chroma namespace.
//
//  Revision 1.19  2004/07/28 02:38:04  edwards
//  Changed {START,END}_CODE("foo") to {START,END}_CODE().
//
//  Revision 1.18  2004/02/11 12:51:34  bjoo
//  Stripped out Read() and Write()
//
//  Revision 1.17  2004/02/03 20:47:24  edwards
//  Small code tweaks.
//
//  Revision 1.16  2003/10/01 03:01:39  edwards
//  Removed extraneous include.
//
//  Revision 1.15  2003/09/29 21:31:36  edwards
//  Tiny cosmetic change.
//
//  Revision 1.14  2003/06/24 03:25:27  edwards
//  Changed from nml to xml.
//
//  Revision 1.13  2003/04/02 22:28:22  edwards
//  Changed proto.h to qdp_util.h
//
//  Revision 1.12  2003/04/01 03:27:26  edwards
//  Added const to sftmom.
//
//  Revision 1.11  2003/04/01 02:38:26  edwards
//  Added doxygen comments.
//
//  Revision 1.10  2003/03/14 21:51:54  flemingg
//  Changes the way in which the nml data is output to match what's done
//  in szin.
//
//  Revision 1.9  2003/03/14 17:16:13  flemingg
//  Variant 1 is now working with SftMom::sft().  In arbitrary units,
//  the relative performance seems to be: V1) 7.5, V2) 10, V3) 100.
//
//  Revision 1.8  2003/03/14 05:14:32  flemingg
//  rewrite of mesons_w.cc to use the new SftMom class.  mesons_w.cc still
//  needs to be cleaned up once the best strategy is resolved.  But for now,
//  the library and test program compiles and runs.
//
//  Revision 1.7  2003/03/06 03:38:35  edwards
//  Added start/end_code.
//
//  Revision 1.6  2003/03/06 02:07:12  flemingg
//  Changed the MomList class to eliminate an unneeded class member.
//
//  Revision 1.5  2003/03/06 00:30:14  flemingg
//  Complete rewrite of lib/meas/hadron/mesons_w.cc, including a simple test
//  program in mainprogs/tests built with 'make check' and various other
//  changes to autoconf/make files to support this rewrite.
//

#include "chromabase.h"
#include "util/ft/sftmom.h"
#include "meas/hadron/mesons2_qcdsf_w.h"
#include "meas/glue/mesfield.h"

namespace Chroma {

  void write(BinaryWriter& bin, const Mesons_corr_QCDSF_t& mes)
  {
    //write( bin , mes.insert_mom );
    write( bin , mes.correlator , mes.correlator.size() );
  }
  void write(BinaryWriter& bin, const Mesons_mom_QCDSF_t& mes)
  {
    write( bin , mes.momentum , mes.momentum.size() );
  }
  void write(BinaryWriter& bin, const Mesons_gamma2_QCDSF_t& mes)
  {
    write( bin , mes.gamma_value , mes.gamma_value.size() );
  }
  void write(BinaryWriter& bin, const MesonsQCDSF_t& mesons)
  {
    write( bin , mesons.gamma_value , mesons.gamma_value.size() );
  }



  //! Meson 2-pt functions
  /* This routine is specific to Wilson fermions!
   *
   * Construct meson propagators and writes in COMPLEX 
   *
   * The two propagators can be identical or different.
   *
   * \param quark_prop_1  first quark propagator ( Read )
   * \param quark_prop_2  second (anti-) quark propagator ( Read )
   * \param t0            timeslice coordinate of the source ( Read )
   * \param phases        object holds list of momenta and Fourier phases ( Read )
   * \param xml           xml file object ( Write )
   * \param xml_group     string used for writing xml data ( Read )
   *
   *        ____
   *        \
   * m(t) =  >  < m(t_source, 0) m(t + t_source, x) >
   *        /
   *        ----
   *          x
   */


  void mesons2qcdsf(const LatticePropagator& quark_prop_1,
		    const LatticePropagator& quark_prop_2,
		    const SftMom& phases,
		    int t0,
		    MesonsQCDSF_t& mesons)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets();

    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    // This variant uses the function SftMom::sft() to do all the work
    // computing the Fourier transform of the meson correlation function
    // inside the class SftMom where all the of the Fourier phases and
    // momenta are stored.  It's primary disadvantage is that it
    // requires more memory because it does all of the Fourier transforms
    // at the same time.

    // Loop over gamma matrix insertions

    mesons.gamma_value.resize(Ns*Ns);

    for (int gamma2=0; gamma2 < (Ns*Ns); ++gamma2)
      {
	mesons.gamma_value[gamma2].gamma_value.resize(Ns*Ns);
	for (int gamma1=0; gamma1 < (Ns*Ns); ++gamma1)
	  {
	    // push(xml_gamma);     // next array element
	    // write(xml_gamma, "gamma_value", gamma2);
	    // write(xml_gamma, "gamma_value", gamma1);

	    // Construct the meson correlation function
	    LatticeComplex corr_fn;
	    corr_fn = trace(adj(anti_quark_prop) * (Gamma(gamma2) *
						    quark_prop_1 * Gamma(gamma1)));

	    multi2d<DComplex> hsum;
	    hsum = phases.sft(corr_fn);

	    // Loop over sink momenta
	    // XMLArrayWriter xml_sink_mom(xml_gamma,phases.numMom());
	    // push(xml_sink_mom, "momenta");

	    mesons.gamma_value[gamma2].gamma_value[gamma1].momentum.resize(phases.numMom());

	    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	      {
		// push(xml_sink_mom);
		// write(xml_sink_mom, "sink_mom_num", sink_mom_num);
		// write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

		mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].correlator.resize(length);
		//mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].insert_mom = phases.numToMom(sink_mom_num);

		//multi1d<DComplex> mesprop(length);
		for (int t=0; t < length; ++t) 
		  {
		    int t_eff = (t - t0 + length) % length;
		    //mesprop[t_eff] = hsum[sink_mom_num][t];
		    mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].correlator[t_eff] = hsum[sink_mom_num][t];
		  }

		//write(xml_sink_mom, "mesprop", mesprop);
		//pop(xml_sink_mom);

	      } // end for(sink_mom_num)
 
	    //pop(xml_sink_mom);
	    //pop(xml_gamma);
	  } // end for(gamma_value)
      }
    //pop(xml_gamma);

    END_CODE();
  }


  void concur2qcdsf(const multi1d<LatticeColorMatrix>& u, 
		    const LatticePropagator& quark_prop_1,
		    const LatticePropagator& quark_prop_2,
		    const SftMom& phases,
		    int t0,
		    MesonsQCDSF_t& mesons)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets();

    LatticePropagator tmp_prop1;
    LatticePropagator tmp_prop2;
    LatticePropagator tmp_prop3;
    LatticeReal psi_sq;
    LatticeReal chi_sq;
    multi1d<LatticeColorMatrix> tfmunu;

    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    // This variant uses the function SftMom::sft() to do all the work
    // computing the Fourier transform of the meson correlation function
    // inside the class SftMom where all the of the Fourier phases and
    // momenta are stored.  It's primary disadvantage is that it
    // requires more memory because it does all of the Fourier transforms
    // at the same time.

    // Loop over gamma matrix insertions

    mesons.gamma_value.resize(Ns*Ns);

    for (int gamma2=0; gamma2 < (Ns*Ns); ++gamma2)
      {
	mesons.gamma_value[gamma2].gamma_value.resize(Ns*Ns);

	int mu, tgamma, tgamma2, offset;

	switch(gamma2){
	case  1:
	case 14:
	  mu = 0;
	  tgamma=gamma2;
	  break;
	case  2:
	case 13:
	  mu = 1;
	  tgamma=gamma2;
	  break;
	case  4:
	case 11:
	  mu = 2;
	  tgamma=gamma2;
	  break;
	case  8:
	case  7:
	  mu = 3;
	  tgamma=gamma2;
	  break;
	case 9:
	  mu = 0;
	  tgamma=14;
	  break;
	case 10:
	  mu = 1;
	  tgamma=13;
	  break;
	case 12:
	  mu = 2;
	  tgamma=11;
	  break;
	case  3:
	  mu = 3;
	  tgamma=7;
	  break;
	default:
	  mu = -1;
	}

	switch(gamma2){
	case 14:
	case 13:
	case 11:
	case 7:
	  tgamma2=15;
	  break;
	default:
	  tgamma2=0;
	}

	switch(gamma2){
	case  1:
	case 14:
	case  2:
	case 13:
	case  4:
	case 11:
	case  8:
	case  7:
	  tmp_prop1 = u[mu] * shift(quark_prop_1, FORWARD, mu);
	  tmp_prop2 = 0.5 * adj(anti_quark_prop) * ( Gamma(tgamma2) * tmp_prop1  + Gamma(tgamma) * tmp_prop1 );

	  tmp_prop1 = u[mu] * shift(anti_quark_prop, FORWARD, mu);
	  
	  tmp_prop2 -= 0.5 * adj(tmp_prop1) * ( Gamma(tgamma2) * quark_prop_1  - Gamma(tgamma) * quark_prop_1 );
	  break;

	case 9:
	case 10:
	case 12:
	case  3:
	  tmp_prop1 = u[mu] * shift(quark_prop_1, FORWARD, mu);
	  tmp_prop2 = 0.5 * adj(anti_quark_prop) * Gamma(tgamma) * tmp_prop1;

	  tmp_prop1 = u[mu] * shift(anti_quark_prop, FORWARD, mu);
	  
	  tmp_prop2 += 0.5 * adj(tmp_prop1) * Gamma(tgamma) * quark_prop_1;
	  break;

	case 5:
	  // Wilson term
	case 6:
	  if (gamma2==5) {
	    tmp_prop2 = 0;
	  }
	  else if (gamma2==6) {
	    tmp_prop2 = 
	      8. * adj(anti_quark_prop) * Gamma(15) * quark_prop_1;
	  }
	  for (int nu=0; nu < Nd ; ++nu) 
	    {	 
	      if (gamma2==5) {
		
		if (nu==0)  tgamma=14;
		else if (nu==1) tgamma=13;
		else if (nu==2) tgamma=11;
		else if (nu==3) tgamma=7;
	      }
	      else if (gamma2==6) tgamma=15;

	      tmp_prop1 = u[nu] * shift(quark_prop_1, FORWARD, nu);
	      if (gamma2==5) {
		tmp_prop2 += 0.5 * adj(anti_quark_prop) * Gamma(tgamma) * tmp_prop1;
	      }
	      else if (gamma2==6) {
		tmp_prop2 -= 0.5 * adj(anti_quark_prop) * Gamma(tgamma) * tmp_prop1;
	      }

	      tmp_prop1 = u[nu] * shift(anti_quark_prop, FORWARD, nu);
	      if (gamma2==5) {
		tmp_prop2 += 0.5 * adj(tmp_prop1) * Gamma(tgamma) * quark_prop_1;
	      }
	      else if (gamma2==6) {
		tmp_prop2 -= 0.5 * adj(tmp_prop1) * Gamma(tgamma) * quark_prop_1;
	      }

	      tmp_prop1 = shift(u[nu], BACKWARD, nu) * quark_prop_1;
	      tmp_prop2 -= 0.5 * adj(shift(anti_quark_prop, BACKWARD, nu)) * Gamma(tgamma) * tmp_prop1;

	      tmp_prop1 = shift(u[nu], BACKWARD, nu) * anti_quark_prop;
	      tmp_prop2 -= 0.5 * adj(tmp_prop1) * Gamma(tgamma) * shift(quark_prop_1, BACKWARD, nu);
	    }

	  break;
	  //Clover term (not yet coded, just scalar for now)
	case 0:
	  mesField(tfmunu, u);

	  offset = 0;
	  tmp_prop2 = 0;
	  // -2 *( CSW * Sum_mu<nu (sigmamunu*Fmunu) / 4 )
	  // sigmamunu = 2*gammamu*gammanu
	  //NOTE: TODO: Remove this extra factor of 2 (i.e. want O_C not 2*O_C)
	  for(int imu=0; imu < Nd-1; ++imu)
	    {
	      int nmu = 1 << imu;
	      for(int inu=imu+1; inu < Nd; ++inu)
		{
		  int nnu = 1 << inu;

		  tmp_prop2 -= adj(anti_quark_prop) * Gamma(15) * Gamma(nmu) * Gamma(nnu) * tfmunu[offset]* quark_prop_1;

		  ++offset;
		}
	    }
	  tmp_prop2 *= 2.65;
	  break;

	  //Pseudoscalar term
	case 15:
	  tmp_prop2 = adj(anti_quark_prop) * Gamma(gamma2) * quark_prop_1;
	  break;
	}

	for (int gamma1=0; gamma1 < (Ns*Ns); ++gamma1)
	  {
	    // push(xml_gamma);     // next array element
	    // write(xml_gamma, "gamma_value", gamma2);
	    // write(xml_gamma, "gamma_value", gamma1);

	    // Construct the meson correlation function
	    LatticeComplex corr_fn;
	    corr_fn = trace( tmp_prop2 * Gamma(gamma1));

	    multi2d<DComplex> hsum;
	    hsum = phases.sft(corr_fn);

	    // Loop over sink momenta
	    // XMLArrayWriter xml_sink_mom(xml_gamma,phases.numMom());
	    // push(xml_sink_mom, "momenta");

	    mesons.gamma_value[gamma2].gamma_value[gamma1].momentum.resize(phases.numMom());

	    for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	      {
		// push(xml_sink_mom);
		// write(xml_sink_mom, "sink_mom_num", sink_mom_num);
		// write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num));

		mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].correlator.resize(length);
		//mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].insert_mom = phases.numToMom(sink_mom_num);

		//multi1d<DComplex> mesprop(length);
		for (int t=0; t < length; ++t) 
		  {
		    int t_eff = (t - t0 + length) % length;
		    //mesprop[t_eff] = hsum[sink_mom_num][t];
		    mesons.gamma_value[gamma2].gamma_value[gamma1].momentum[sink_mom_num].correlator[t_eff] = hsum[sink_mom_num][t];
		  }

		//write(xml_sink_mom, "mesprop", mesprop);
		//pop(xml_sink_mom);

	      } // end for(sink_mom_num)
 
	    //pop(xml_sink_mom);
	    //pop(xml_gamma);
	  } // end for(gamma_value)
      }
    //pop(xml_gamma);

    END_CODE();
  }


  void mesons2qcdsfsmall(const LatticePropagator& quark_prop_1,
			 const LatticePropagator& quark_prop_2,
			 const SftMom& phases,
			 int t0,
			 Mesons_gamma2_QCDSF_t& mesons)
  {
    START_CODE();
    
    int ngamma=20;

    // list of gammas to loop over
    multi1d<int>     gammalist1(ngamma);
    multi1d<int>     gammalist2(ngamma);

    gammalist1[0]=15; gammalist2[0]=15;
    gammalist1[1]=15; gammalist2[1]=7;
    gammalist1[2]=7;  gammalist2[2]=15;
    gammalist1[3]=7;  gammalist2[3]=7;
    gammalist1[4]=1;  gammalist2[4]=1;
    gammalist1[5]=1;  gammalist2[5]=2;
    gammalist1[6]=1;  gammalist2[6]=4;
    gammalist1[7]=1;  gammalist2[7]=8;
    gammalist1[8]=2;  gammalist2[8]=1;
    gammalist1[9]=2;  gammalist2[9]=2;
    gammalist1[10]=2; gammalist2[10]=4;
    gammalist1[11]=2; gammalist2[11]=8;
    gammalist1[12]=4; gammalist2[12]=1;
    gammalist1[13]=4; gammalist2[13]=2;
    gammalist1[14]=4; gammalist2[14]=4;
    gammalist1[15]=4; gammalist2[15]=8;
    gammalist1[16]=8; gammalist2[16]=1;
    gammalist1[17]=8; gammalist2[17]=2;
    gammalist1[18]=8; gammalist2[18]=4;
    gammalist1[19]=8; gammalist2[19]=8;

    // Length of lattice in decay direction
    int length = phases.numSubsets();

    // Construct the anti-quark propagator from quark_prop_2
    int G5 = Ns*Ns-1;
    LatticePropagator anti_quark_prop =  Gamma(G5) * quark_prop_2 * Gamma(G5);

    // This variant uses the function SftMom::sft() to do all the work
    // computing the Fourier transform of the meson correlation function
    // inside the class SftMom where all the of the Fourier phases and
    // momenta are stored.  It's primary disadvantage is that it
    // requires more memory because it does all of the Fourier transforms
    // at the same time.

    // Loop over gamma matrix insertions

    mesons.gamma_value.resize(ngamma);

    for (int igamma=0; igamma < ngamma; ++igamma)
      {
	int gamma1=gammalist1[igamma];
	int gamma2=gammalist2[igamma];

	// Construct the meson correlation function
	LatticeComplex corr_fn;
	corr_fn = trace(adj(anti_quark_prop) * (Gamma(gamma2) *
						quark_prop_1 * Gamma(gamma1)));

	multi2d<DComplex> hsum;
	hsum = phases.sft(corr_fn);

	// Loop over sink momenta

	mesons.gamma_value[igamma].momentum.resize(phases.numMom());

	for (int sink_mom_num=0; sink_mom_num < phases.numMom(); ++sink_mom_num) 
	  {
	    
	    mesons.gamma_value[igamma].momentum[sink_mom_num].correlator.resize(length);

	    for (int t=0; t < length; ++t) 
	      {
		int t_eff = (t - t0 + length) % length;
		mesons.gamma_value[igamma].momentum[sink_mom_num].correlator[t_eff] = hsum[sink_mom_num][t];
	      }

	  } // end for(sink_mom_num)
 
      } // end for(gamma_value)

    END_CODE();
  }


}  // end namespace Chroma
