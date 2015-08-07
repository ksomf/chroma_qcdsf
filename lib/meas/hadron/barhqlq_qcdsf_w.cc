// $Id: barhqlq_w.cc,v 3.5 2009-03-19 17:17:20 mcneile Exp $
/*! \file
 *  \brief Heavy-light baryon 2-pt functions
 */


#include "meas/hadron/barhqlq_w.h"
#include "meas/hadron/barhqlq_qcdsf_w.h"
#include "meas/hadron/barspinmat_w.h"

using namespace Chroma;

namespace Chroma
{

  void write(BinaryWriter& bin, const Baryons_corr_QCDSF_t& mes)
  {
    //write( bin , mes.insert_mom );
    write( bin , mes.correlator , mes.correlator.size() );
  }
  void write(BinaryWriter& bin, const Baryons_N_QCDSF_t& mes)
  {
    write( bin , mes.momentum , mes.momentum.size() );
  }
  void write(BinaryWriter& bin, const BaryonsQCDSF_t& baryons)
  {
    write( bin , baryons.barnum , baryons.barnum.size() );
  }


  void barhqlq_qcdsf_lime(const LatticePropagator& propagator_1,
			  const LatticePropagator& propagator_2,
			  const LatticePropagator& propagator_3,
			  const bool & haveThird,
			  const SftMom& phases,
			  int t0, int bc_spec, bool time_rev, bool fwdbwd_average,
			  BaryonsQCDSF_t& bar,
			  BaryonsQCDSF_t& bar_trev)
  {
    START_CODE();

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    multi3d<DComplex> bardisp1;
    multi3d<DComplex> bardisp2;

    // Forward
    barhqlq_qcdsf(propagator_1, propagator_2, propagator_3, haveThird, phases, bardisp1);

    // Possibly add in a time-reversed contribution
    bool time_revP = (bc_spec*bc_spec == 1) ? time_rev : false;

    if (time_revP) {
      if (haveThird) {
	/* Time-charge reverse the quark propagators */
	/* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
	LatticePropagator q1_tmp = - (Gamma(7) * propagator_1 * Gamma(7));
	LatticePropagator q2_tmp = - (Gamma(7) * propagator_2 * Gamma(7));
	LatticePropagator q3_tmp = - (Gamma(7) * propagator_3 * Gamma(7));

	barhqlq_qcdsf(q1_tmp, q2_tmp, q3_tmp, haveThird, phases, bardisp2);
      } else {
	/* Time-charge reverse the quark propagators */
	/* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
	LatticePropagator q1_tmp = - (Gamma(7) * propagator_1 * Gamma(7));
	LatticePropagator q2_tmp = - (Gamma(7) * propagator_2 * Gamma(7));

	barhqlq_qcdsf(q1_tmp, q2_tmp, q2_tmp, haveThird, phases, bardisp2);
      }
    }

    int num_baryons = bardisp1.size3();
    int num_mom = bardisp1.size2();
    int length  = bardisp1.size1();

    // Loop over baryons
    // XMLArrayWriter xml_bar(xml,num_baryons);
    // push(xml_bar, xml_group);

    bar.barnum.resize(num_baryons);
    bar_trev.barnum.resize(num_baryons);
    for(int baryons = 0; baryons < num_baryons; ++baryons)
      {
	// push(xml_bar);     // next array element
	// write(xml_bar, "baryon_num", baryons);

	// // Loop over sink momenta
	// XMLArrayWriter xml_sink_mom(xml_bar,num_mom);
	// push(xml_sink_mom, "momenta");

	bar.barnum[baryons].momentum.resize(num_mom);
	bar_trev.barnum[baryons].momentum.resize(num_mom);
	for(int sink_mom_num = 0; sink_mom_num < num_mom; ++sink_mom_num)
	  {
	    // push(xml_sink_mom);
	    // write(xml_sink_mom, "sink_mom_num", sink_mom_num) ;
	    // write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num)) ;

	    multi1d<Complex> barprop(length);

	    /* forward */
	    for(int t = 0; t < length; ++t)
	      {
		int t_eff = (t - t0 + length) % length;

		if ( bc_spec < 0 && (t_eff+t0) >= length)
		  barprop[t_eff] = -bardisp1[baryons][sink_mom_num][t];
		else
		  barprop[t_eff] =  bardisp1[baryons][sink_mom_num][t];
	      }

	    multi1d<Complex> barprop_trev(length);

	    if (time_revP)
	      {
		/* backward */
		for(int t = 0; t < length; ++t)
		  {
		    int t_eff = (length - t + t0) % length;

		    if ( bc_spec < 0 && (t_eff-t0) > 0)
		      {
			barprop_trev[t_eff] = -bardisp2[baryons][sink_mom_num][t];
		      }
		    else
		      {
			barprop_trev[t_eff] = bardisp2[baryons][sink_mom_num][t];
		      }
		  }

		if(fwdbwd_average){
		  barprop += barprop_trev;
		  barprop *= 0.5;
		}
	      }

	    //bardisp1[baryons][sink_mom_num] = barprop;
	    bar.barnum[baryons].momentum[sink_mom_num].correlator.resize(length);
	    bar_trev.barnum[baryons].momentum[sink_mom_num].correlator.resize(length);

	    bar.barnum[baryons].momentum[sink_mom_num].correlator = barprop;
	    bar_trev.barnum[baryons].momentum[sink_mom_num].correlator = barprop_trev;

	    // write(xml_sink_mom, "barprop", barprop);
	    // pop(xml_sink_mom);
	  } // end for(sink_mom_num)

	// pop(xml_sink_mom);
	// pop(xml_bar);
      } // end for(gamma_value)

    // pop(xml_bar);

    END_CODE();
  }



  // Original Chroma code (timerev behaviour)
	    // multi1d<Complex> barprop(length);

	    // /* forward */
	    // for(int t = 0; t < length; ++t)
	    //   {
	    // 	int t_eff = (t - t0 + length) % length;

	    // 	if ( bc_spec < 0 && (t_eff+t0) >= length)
	    // 	  barprop[t_eff] = -bardisp1[baryons][sink_mom_num][t];
	    // 	else
	    // 	  barprop[t_eff] =  bardisp1[baryons][sink_mom_num][t];
	    //   }

	    // if (time_revP)
	    //   {
	    // 	/* backward */
	    // 	for(int t = 0; t < length; ++t)
	    // 	  {
	    // 	    int t_eff = (length - t + t0) % length;

	    // 	    if ( bc_spec < 0 && (t_eff-t0) > 0)
	    // 	      {
	    // 		barprop[t_eff] -= bardisp2[baryons][sink_mom_num][t];
	    // 		barprop[t_eff] *= 0.5;
	    // 	      }
	    // 	    else
	    // 	      {
	    // 		barprop[t_eff] += bardisp2[baryons][sink_mom_num][t];
	    // 		barprop[t_eff] *= 0.5;
	    // 	      }
	    // 	  }
	    //   }




  void barhqlq_qcdsf_xml(const LatticePropagator& propagator_1,
			 const LatticePropagator& propagator_2,
			 const LatticePropagator& propagator_3,
			 const bool & haveThird,
			 const SftMom& phases,
			 int t0, int bc_spec, bool time_rev, bool fwdbwd_average,
			 XMLWriter& xml,
			 const std::string& xml_group)
  {
    START_CODE();

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    multi3d<DComplex> bardisp1;
    multi3d<DComplex> bardisp2;

    // Forward
    barhqlq_qcdsf(propagator_1, propagator_2,propagator_3,haveThird, phases, bardisp1);

    // Possibly add in a time-reversed contribution
    bool time_revP = (bc_spec*bc_spec == 1) ? time_rev : false;

    if (time_revP) {
      if (haveThird) {
	/* Time-charge reverse the quark propagators */
	/* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
	LatticePropagator q1_tmp = - (Gamma(7) * propagator_1 * Gamma(7));
	LatticePropagator q2_tmp = - (Gamma(7) * propagator_2 * Gamma(7));
	LatticePropagator q3_tmp = - (Gamma(7) * propagator_3 * Gamma(7));

	barhqlq_qcdsf(q1_tmp, q2_tmp, q3_tmp, haveThird, phases, bardisp2);
      } else {
	/* Time-charge reverse the quark propagators */
	/* S_{CT} = gamma_5 gamma_4 = gamma_1 gamma_2 gamma_3 = Gamma(7) */
	LatticePropagator q1_tmp = - (Gamma(7) * propagator_1 * Gamma(7));
	LatticePropagator q2_tmp = - (Gamma(7) * propagator_2 * Gamma(7));

	barhqlq_qcdsf(q1_tmp, q2_tmp, q2_tmp, false, phases, bardisp2);
      }
    }


    int num_baryons = bardisp1.size3();
    int num_mom = bardisp1.size2();
    int length  = bardisp1.size1();

    // Loop over baryons
    XMLArrayWriter xml_bar(xml,num_baryons);
    push(xml_bar, xml_group);

    for(int baryons = 0; baryons < num_baryons; ++baryons)
    {
      push(xml_bar);     // next array element
      write(xml_bar, "baryon_num", baryons);

      // Loop over sink momenta
      XMLArrayWriter xml_sink_mom(xml_bar,num_mom);
      push(xml_sink_mom, "momenta");

      for(int sink_mom_num = 0; sink_mom_num < num_mom; ++sink_mom_num)
      {
	push(xml_sink_mom);
	write(xml_sink_mom, "sink_mom_num", sink_mom_num) ;
	write(xml_sink_mom, "sink_mom", phases.numToMom(sink_mom_num)) ;

	multi1d<Complex> barprop(length);

	/* forward */
	for(int t = 0; t < length; ++t)
	{
	  int t_eff = (t - t0 + length) % length;

	  if ( bc_spec < 0 && (t_eff+t0) >= length)
	    barprop[t_eff] = -bardisp1[baryons][sink_mom_num][t];
	  else
	    barprop[t_eff] =  bardisp1[baryons][sink_mom_num][t];
	}

	multi1d<Complex> barprop_trev(length);

	if (time_revP)
	{
	  /* backward */
	  for(int t = 0; t < length; ++t)
	  {
	    int t_eff = (length - t + t0) % length;

	    if ( bc_spec < 0 && (t_eff-t0) > 0)
	    {
	      barprop_trev[t_eff] = -bardisp2[baryons][sink_mom_num][t];
	    }
	    else
	    {
	      barprop_trev[t_eff] = bardisp2[baryons][sink_mom_num][t];
	    }
	  }

	  if(fwdbwd_average){
	    barprop += barprop_trev;
	    barprop *= 0.5;
	  }
	}

	write(xml_sink_mom, "barprop", barprop);
	if(time_revP && !fwdbwd_average){
	  write(xml_sink_mom, "barprop_trev", barprop_trev);
	}

	pop(xml_sink_mom);
      } // end for(sink_mom_num)

      pop(xml_sink_mom);
      pop(xml_bar);
    } // end for(gamma_value)

    pop(xml_bar);

    END_CODE();
  }


  void barhqlq_qcdsf(const LatticePropagator& quark_propagator_1,
		     const LatticePropagator& quark_propagator_2,
		     const LatticePropagator& quark_propagator_3,
		     const bool & haveThird,
		     const SftMom& phases,
		     multi3d<DComplex>& barprop)
  {
    START_CODE();

    // Length of lattice in decay direction
    int length = phases.numSubsets() ;

    if ( Ns != 4 || Nc != 3 )		/* Code is specific to Ns=4 and Nc=3. */
      return;

    // Setup the return stuff
    int num_baryons_help;
    if (!haveThird) num_baryons_help= 26;
    else num_baryons_help=36;

    const int num_baryons=num_baryons_help;
    int num_mom = phases.numMom();
    barprop.resize(num_baryons,num_mom,length);


    // T_mixed = (1 + \Sigma_3)*(1 + gamma_4) / 2
    //         = (1 + Gamma(8) - i G(3) - i G(11)) / 2
    SpinMatrix T_mixed = BaryonSpinMats::Tmixed();

    // T_mixed = (1 - \Sigma_3)*(1 + gamma_4) / 2
    //         = (1 + Gamma(8) + i G(3) + i G(11)) / 2
    SpinMatrix T_mixedminus = BaryonSpinMats::Tmixedminus();

    // T_unpol = (1/2)(1 + gamma_4)
    SpinMatrix T_unpol = BaryonSpinMats::Tunpol();

    // T_unpolg5 = (1/2)(1 + gamma_4) * g5
    SpinMatrix T_unpolg5 = BaryonSpinMats::Tunpolg5();

    //! T = \Sigma_1 (1 + gamma_4) / 2
    SpinMatrix Tpolx = BaryonSpinMats::Tpolx();

    //! T = \Sigma_2 (1 + gamma_4) / 2
    SpinMatrix Tpoly = BaryonSpinMats::Tpoly();

    //! T = \Sigma_3 (1 + gamma_4) / 2
    SpinMatrix Tpol = BaryonSpinMats::Tpol();

    // C gamma_5 = Gamma(5)
    SpinMatrix Cg5 = BaryonSpinMats::Cg5();

    // C gamma_5 gamma_4 = - Gamma(13)
    SpinMatrix Cg5g4 = BaryonSpinMats::Cg5g4();

    // C g_5 NR = (1/2)*C gamma_5 * ( 1 + g_4 )
    SpinMatrix Cg5NR = BaryonSpinMats::Cg5NR();

    LatticeComplex b_prop;

    // Loop over baryons
    for(int baryons = 0; baryons < num_baryons; ++baryons)
    {
      LatticePropagator di_quark;

      switch (baryons)
      {
      case 0:
	// Sigma^+_1 (or proton); use also for Lambda_1!
	// |S_1, s_z=1/2> = (s C gamma_5 u) "u_up"
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_mixed, Cg5);
	break;

      case 1:
	// Lambda_1
	// |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
	//                  + (u C gamma_5 s) "d_up" , see comments at top
	// C gamma_5 = Gamma(5)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions_3prop::lambda8_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_2,
						  T_mixed, Cg5);
	break;

      case 2:
	// Sigma^{*+}_1
	// |S*_1, s_z=3/2> = 2*(s C gamma_- u) "u_up" + (u C gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixed, BaryonSpinMats::Cgm());
	break;

      case 3:
	// Sigma^+_2; use also for Lambda_2!
	// |S_2, s_z=1/2> = (s C gamma_4 gamma_5 u) "u_up"
	// C gamma_5 gamma_4 = - Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_mixed, Cg5g4);
	break;

      case 4:
	// Lambda_2
	// |L_2, s_z=1/2> = 2*(u C gamma_4 gamma_5 d) "s_up"
	//                  + (s C gamma_4 gamma_5 d) "u_up"
	//                  + (u C gamma_4 gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions_3prop::lambda8_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_2, T_mixed, Cg5g4);
	break;

      case 5:
	// Sigma^{*+}_2
	// |S*_2, s_z=3/2> = 2*(s C gamma_4 gamma_- u) "u_up" + (u C gamma_4 gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixed, BaryonSpinMats::Cg4m());
	break;

      case 6:
	// Sigma^+_3; use also for Lambda_3!
	// |S_3, s_z=1/2> = (s C (1/2)(1 + gamma_4) gamma_5 u) "u_up"
	// C gamma_5 gamma_4 = - Gamma(13)
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_mixed, Cg5NR);
	break;

      case 7:
	// Lambda_3
	// |L_3, s_z=1/2> = 2*(u C (1/2)(1 + gamma_4) gamma_5 d) "s_up"
	//                  + (s C (1/2)(1 + gamma_4) gamma_5 d) "u_up"
	//                  + (u C (1/2)(1 + gamma_4) gamma_5 s) "d_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions_3prop::lambda8_2pt(quark_propagator_1, quark_propagator_2,quark_propagator_2,
						  T_mixed, Cg5NR);
	break;

      case 8:
	// Sigma^{*+}_3
	// |S*_3, s_z=3/2> = 2*(s C (1/2)(1+gamma_4) gamma_- u) "u_up"
	//                   + (u C (1/2)(1+gamma_4) gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	// Arrgh, goofy CgmNR normalization again from szin code.
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixed, BaryonSpinMats::CgmNR());

	// Agghh, we have a goofy factor of 4 normalization factor here. The
	// ancient szin way didn't care about norms, so it happily made it
	// 4 times too big. There is a missing 0.5 in the NR normalization
	// in the old szin code.
	// So, we compensate to keep the same normalization
	b_prop *= 4.0;
	break;


      case 9:
	// Sigma^+_4 -- but unpolarised
	// |S_4, s_z=1/2> = (s C gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_unpol, Cg5);
	break;

      case 10:
	// Sigma^+_5
	// |S_5, s_z=1/2> = (s C gamma_4 gamma_5 u) "u_up", see comments at top
	// C gamma_5 gamma_4 = - Gamma(13)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_unpol, Cg5g4);
	break;

      case 11:
	// Sigma^+_6
	// |S_6, s_z=1/2> = (s C (1/2)(1 + gamma_4) gamma_5 u) "u_up", see comments at top
	// C gamma_5 = Gamma(5)
	// Unpolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 T_unpol, Cg5NR);
	break;

      case 12:
	// Lambda_4 : naive Lambda interpolating field
	// |L_4 > = (d C gamma_5 u) s
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::lambdaNaive2pt(quark_propagator_1, quark_propagator_2,
						       T_unpol, Cg5);
	break;

      case 13:
	// Xi_1
	// |X_1 > = (s C gamma_5 u) s
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions::xi2pt(quark_propagator_1, quark_propagator_2,
					      T_unpol, Cg5);
	break;

      case 14:
	// Lambda_5 : naive Lambda interpolating field
	// |L_5 > = (d C gamma_5 u) "s_up"
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::lambdaNaive2pt(quark_propagator_1, quark_propagator_2,
						       T_unpol, Cg5);
	break;

      case 15:
	// Xi_2
	// |X_2 > = (s C gamma_5 u) "s_up"
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::xi2pt(quark_propagator_1, quark_propagator_2,
					      T_mixed, Cg5);
	break;

      case 16:
	// Proton_negpar_3; use also for Lambda_negpar_3!
	// |P_7, s_z=1/2> = (d C gamma_5 (1/2)(1 - g_4) u) "u_up", see comments at top
	// C g_5 NR negpar = (1/2)*C gamma_5 * ( 1 - g_4 )
	// T = (1 + \Sigma_3)*(1 - gamma_4) / 2
	//   = (1 - Gamma(8) + i G(3) - i G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
						 BaryonSpinMats::TmixedNegPar(),
						 BaryonSpinMats::Cg5NRnegPar());
	break;
      case 17:
	// Proton_pol_1;
	// C g_5 = gamma_5
	// T = \Sigma_1*(1 + gamma_4) / 2
	//   = -i(G(6) + G(14)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, Tpolx, Cg5);
	break;
      case 18:
	// Proton_pol_2;
	// C g_5 = gamma_5
	// T = \Sigma_2*(1 + gamma_4) / 2
	//   = -i(G(5) + G(13)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, Tpoly, Cg5);
	break;
      case 19:
	// Proton_pol_3;
	// C g_5 = gamma_5
	// T = \Sigma_3*(1 + gamma_4) / 2
	//   = -i(G(3) + G(11)) / 2
	b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2, Tpol, Cg5);
	break;

      case 20:
	// Sigma^+_1 (or proton); use also for Lambda_1!
        // |S_1, s_z=1/2> = (s C gamma_5 u) "u_up"
        // C gamma_5 = Gamma(5)
        // Unpolarized but with g5:
        // T_unpolg5 = T = gamma_5*(1 + gamma_4) / 2
        //             = (Gamma(15) + Gamma(7)) / 2
        b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
                                                 T_unpolg5, Cg5);
	break;

      case 21:
	// Sigma^+_2; use also for Lambda_2!
        // This used to be case 0 but my ugly hack above means I have moved it here
        // |S_2, s_z=-1/2> = (s C gamma_5 u) "u_down"
        // Polarized:
        // T_mixedminus = T = (1 - \Sigma_3)*(1 + gamma_4) / 2
        //                  = (1 + Gamma(8) + i G(3) + i G(11)) / 2
        b_prop = Baryon2PtContractions::sigma2pt(quark_propagator_1, quark_propagator_2,
                                                 T_mixedminus, Cg5);
	break;


      case 22:
	// Lambda_1
	// |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_down" + (s C gamma_5 d) "u_down"
	//                  + (u C gamma_5 s) "d_down" , see comments at top
	// C gamma_5 = Gamma(5)
	// Polarized:
        // T_mixedminus = T = (1 - \Sigma_3)*(1 + gamma_4) / 2
        //                  = (1 + Gamma(8) + i G(3) + i G(11)) / 2
	b_prop = Baryon2PtContractions_3prop::lambda8_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_2,
						  T_mixedminus, Cg5);
	break;

      case 23:
	// Sigma^{*+}_1
	// |S*_1, s_z=1/2> = 2*(s C gamma_- u) "u_down" + (u C gamma_- u) "s_down"
	// Polarized:
        // T_mixedminus = T = (1 - \Sigma_3)*(1 + gamma_4) / 2
        //                  = (1 + Gamma(8) + i G(3) + i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixedminus, BaryonSpinMats::Cgm());
	break;

      case 24:
	// Sigma^{*+}_1
	// |S*_1, s_z=-1/2> = 2*(s C gamma_- u) "u_down" + (u C gamma_- u) "s_down"
	// Polarized:
        // T_mixedminus = T = (1 - \Sigma_3)*(1 + gamma_4) / 2
        //                  = (1 + Gamma(8) + i G(3) + i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixed, BaryonSpinMats::Cgp());
	break;

      case 25:
	// Sigma^{*+}_1
	// |S*_1, s_z=-3/2> = 2*(s C gamma_- u) "u_down" + (u C gamma_- u) "s_down"
	// Polarized:
        // T_mixedminus = T = (1 - \Sigma_3)*(1 + gamma_4) / 2
        //                  = (1 + Gamma(8) + i G(3) + i G(11)) / 2
	b_prop = Baryon2PtContractions::sigmast2pt(quark_propagator_1, quark_propagator_2,
						   T_mixedminus, BaryonSpinMats::Cgp());
	break;

     case 26:
	// Sigma_0;
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)

	b_prop = Baryon2PtContractions_3prop::sigma0_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3,
						  T_unpol, Cg5);


	break;

	case 27:
	// Lambda_octet;
	// Lambda_1
	// |L_1, s_z=1/2> = 2*(u C gamma_5 d) "s_up" + (s C gamma_5 d) "u_up"
	//                  + (u C gamma_5 s) "d_up" , see comments at top
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions_3prop::lambda8_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3,
						  T_unpol, Cg5);
	break;

	case 28:
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions_3prop::lambda8_to_sigma0_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3,
						  T_unpol, Cg5);

	break;

	case 29:
	// C gamma_5 = Gamma(5)
	// UnPolarized:
	// T_unpol = T = (1/2)(1 + gamma_4)
	b_prop = Baryon2PtContractions_3prop::sigma0_lambda8_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3,
						  T_unpol, Cg5);
	 break;

	case 30:
	// Sigma^{*+}_3
	// |S*_3, s_z=3/2> = 2*(s C (1/2)(1+gamma_4) gamma_- u) "u_up"
	//                   + (u C (1/2)(1+gamma_4) gamma_- u) "s_up"
	// Polarized:
	// T_mixed = T = (1 + \Sigma_3)*(1 + gamma_4) / 2
	//             = (1 + Gamma(8) - i G(3) - i G(11)) / 2


    //! C gamma_- = Cgm = (C gamma_-)^T
//     SpinMatrix Cgm()
//     {
//       SpinMatrix g_one = 1.0;
//       SpinMatrix g_tmp1 = 0.5 * (Gamma(2) * g_one  +  timesI(Gamma(1) * g_one));
//       return SpinMatrix(Gamma(10) * g_tmp1);
//     }
	b_prop = Baryon2PtContractions_3prop::sigma_star_2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3,
						  T_mixed, BaryonSpinMats::Cgm());
	 break;

 //Now do sequential source with fixed operator rather than fixed sink
              
     case 31:
      //p->n
              
    b_prop = Baryon2PtContractions_3prop::p_to_n2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3, T_unpol, Cg5);
     
     case 32:
     //p->n
              
    b_prop = Baryon2PtContractions_3prop::p_to_n2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3, Tpolx, Cg5);
              
              
     case 33:
     //p->n
              
     b_prop = Baryon2PtContractions_3prop::p_to_n2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3, Tpoly, Cg5);
              
              
     case 34:
     //p->n
              
     b_prop = Baryon2PtContractions_3prop::p_to_n2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3, Tpol, Cg5);
              
              
     case 35:
     //p->p, u quark unpol
              
      b_prop = Baryon2PtContractions_3prop::p_to_p_u2pt(quark_propagator_1, quark_propagator_2, quark_propagator_3, T_unpol, Cg5);
              
     break;
          
              
              
              
              
              
              
      default:
	QDP_error_exit("Unknown baryon", baryons);
      }

      // Project onto zero and if desired non-zero momentum
      multi2d<DComplex> hsum;
      hsum = phases.sft(b_prop);

      for(int sink_mom_num=0; sink_mom_num < num_mom; ++sink_mom_num)
	for(int t = 0; t < length; ++t)
	{
	  // NOTE: there is NO  1/2  multiplying hsum
	  barprop[baryons][sink_mom_num][t] = hsum[sink_mom_num][t];
	}

    } // end loop over baryons

    END_CODE();
  }



}
