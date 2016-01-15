/*! \file npr_vertex_nonlocal_w.cc
* \author Stanislav Kazmin (University Leipzig)
* \date 2016-01-13 13:49:01
 * \brief NPR vertex calculations for the non local axial current.
 */
#include <iomanip>      // std::setprecision

#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_nonlocal_w.h"
#include "meas/hadron/npr_vertex_w.h"

namespace Chroma
{
	
	void BkwdFrwdNonlocal(const LatticePropagator& F, const multi1d<LatticeColorMatrix>& U, QDPFileWriter& qio_file, int& GBB_NLinkPatterns, const multi1d<int>& LinkDirs)
	// TODO (S. Kazmin): linkDirs is not used, do we need it? write it because npr routine need it
	{
		// reset time measurement
		StopWatch TotalTime;
		TotalTime.reset();
		TotalTime.start();
		LatticePropagator B = Gamma(15) * F * Gamma(15); // (S. Kazmin): backwards propagator NOT adj.
		LatticePropagator div_A = Complex(0); //// (S. Kazmin): divergence of A operator
		DPropagator div_A_mean; //// (S. Kazmin): divergence of A operator
		// calculate one link operator in each direction
		for(int mu = 0; mu < Nd; ++mu) // go though all directions // TODO change to all indices with Ns as maximum
		{
			int gamma = 1 << mu; // get gamma value from bits-shift of 1.
			gamma = 15 - gamma; // get the gamma index of the according γ_μ γ_5 operator
			XMLBufferWriter record_xml;
			push(record_xml, "Vertex");
			QDPIO::cout << __func__ << ": LinkDirs = " << LinkDirs << "  gamma = " << gamma << std::endl;
			write(record_xml, "linkDirs", LinkDirs);   // link pattern
			write(record_xml, "gamma", gamma);
			// counts number of link patterns
			GBB_NLinkPatterns++;
			DPropagator Amu_x_mean;
			//// (S. Kazmin): Amu_x is forward directed and is placed at x
			//// (S. Kazmin): A_μ(x) = 1/2 *( B(x) γ_μ γ_5 U_μ(x) F(x+μ) + B(x+μ) adj(U_μ(x)) γ_μ γ_5  F(x))
			//// (S. Kazmin): the 1/2 factor is shifted to later calculations to reduce thew number of used factor*matrix operations 
			LatticePropagator Amu_x = 
			adj(B) * Gamma(gamma) * (U[mu] * shift(F, FORWARD, mu)) 
			+ adj(U[mu] * shift(B, FORWARD, mu)) * Gamma(gamma)  * F;
			//// (S. Kazmin): Amu_xMinusMu is forward directed and is placed at x-μ
			//// (S. Kazmin): A_μ(x-μ) = 1/2 *( B(x-μ) γ_μ γ_5 U_μ(x-μ) F(x) + B(x) adj(U_μ(x-μ)) γ_μ γ_5  F(x-μ))
			LatticePropagator Amu_xMinusMu = 
			adj(shift(B, BACKWARD, mu)) * Gamma(gamma) * (shift(U[mu], BACKWARD, mu) * F) 
			+ adj( shift(U[mu], BACKWARD, mu) * B ) * Gamma(gamma)  * shift(F, BACKWARD, mu);
			// The site's worth of data of interest
			// sum is over the volume at each site
			Amu_x_mean = 0.5 * sum(Amu_x) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at all sites
			
// 					// TODO(S. Kazmin): TEST the propagator entries. Must be removed before big runs.
// 		// print first entries of div_A
// 		std::cout << std::setprecision(2);
// 		std::cout << "div_A(0,0) = \n";
// 		for (int ind = 0; ind < Layout::vol(); ++ ind)
// 		{
// 			std::cout << Amu_x.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
// 		}
// 		std::cout << sum(Amu_x).elem().elem(0,0).elem(0,0).real() << "\n";
// 			//// END TEST
			// append to xml file
			pop(record_xml);
			write(qio_file, record_xml, Amu_x_mean);
			// calculate divergence of A operator
			// LatticePropagator tmp = Amu_x - Amu_xMinusMu;
			// LatticePropagator tmp2 = div_A;
			if (gamma == 7 or gamma == 13) //// (S. Kazmin): because of Gamma convention in chroma
			{
				div_A -= (Amu_x - Amu_xMinusMu); 
				std::cout << "sign change \n";
			}else
			{
				div_A += (Amu_x - Amu_xMinusMu); 
			}
			//div_A = div_A - Amu_xMinusMu;
					std::cout << "div_A(0,0) = \n";
		for (int ind = 0; ind < Layout::vol(); ++ ind)
		{
			std::cout << "A+ = "<< Amu_x.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
			std::cout << "A- = "<< Amu_xMinusMu.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
			std::cout << "div_A = " << div_A.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
			// std::cout << "tmp = " << tmp.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
			// std::cout << "tmp2 = " << tmp2.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
		}
		}
		// calculate divergence of A operator
		div_A_mean = 0.5 * sum(div_A) / Double(Layout::vol()); // normalize by the volume -> mean value of the prop at all sites
		// write divergence of A to XML file as gamma 0 vertex

		// TODO(S. Kazmin): TEST the propagator entries. Must be removed before big runs.
		// print first entries of div_A
		std::cout << std::setprecision(16);
		std::cout << "div_A(0,0) = \n";
		double testSum = 0.;
		for (int ind = 0; ind < Layout::vol(); ++ ind)
		{
						//std::cout << "testSum: " << testSum << " \n";
			testSum += div_A.elem(ind).elem(0, 0).elem(0, 0).real();

			std::cout << div_A.elem(ind).elem(0, 0).elem(0, 0).real() << " \n";
		}
					std::cout << "testSum: " << testSum << " \n";
		LatticePropagator testProp = div_A;
// 		div_A = Real(1);
		std::cout << "Σ(div_A(0,0)) = \n" << testSum << "\n";
		std::cout << div_A_mean.elem().elem(0,0).elem(0,0).real() << "\n";
			for (int ind = 0; ind < Layout::vol(); ++ ind)
			{
				std::cout << std::setprecision(2);
				std::cout << "place = " << ind  << "\n";
				// div_a OUTPUT
				std::cout << "div_A = \n";
				for(int i=0; i<12; i++)
				{
					for(int j=0; j<12; j++)
					{
						std::cout << "(";
						std::cout << div_A.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).real();
						std::cout << " ";
						std::cout << div_A.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).imag();
						std::cout << ")";
					}
					std::cout << "\n";
				}
			}
			//div_A_mean
							std::cout << "div_A_mean = \n";
				for(int i=0; i<12; i++)
				{
					for(int j=0; j<12; j++)
					{
						std::cout << "(";
						std::cout << div_A_mean.elem().elem(i % 4, j % 4).elem(i / 4, j / 4).real();
						std::cout << " ";
						std::cout << div_A_mean.elem().elem(i % 4, j % 4).elem(i / 4, j / 4).imag();
						std::cout << ")";
					}
					std::cout << "\n";
				}
			//// END TEST
					XMLBufferWriter record_xml;
		push(record_xml, "Vertex");
		QDPIO::cout << __func__ << ": LinkDirs = " << LinkDirs << "  gamma = " << 0 << std::endl;
		QDPIO::cout << __func__ << ": write div(A) to Vertex with gamma = " << 0 << std::endl;
		write(record_xml, "linkDirs", LinkDirs);   // link pattern
		write(record_xml, "annotation", "div(A)");   // link pattern
		write(record_xml, "gamma", 0);
		pop(record_xml);
		write(qio_file, record_xml, div_A_mean);
		// print elapsed time
		TotalTime.stop();
		QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
	
	//! NPR vertices
	void NprVertexNonlocal(const LatticePropagator& F,const multi1d< LatticeColorMatrix >& U,const unsigned short int MaxNLinks, const BBLinkPattern LinkPattern, QDPFileWriter& qio_file)
	{
		// start time measurement
		//	StopWatch TotalTime;
		//	TotalTime.reset();
		//	TotalTime.start();
		StopWatch Timer;
		Timer.reset();
		Timer.start();
		int GBB_NLinkPatterns;
		// calculate building blocks
		QDPIO::cout << __func__ << ": start BkwdFrwdNonlocal" << std::endl;
		const int NLinks = 0;
		multi1d<int> LinkDirs(0);
		//// (S. Kazmin): Backward propagator is calculated in BkwdFrwdNonlocal function
		BkwdFrwdNonlocal(F, U, qio_file, GBB_NLinkPatterns, LinkDirs);
		Timer.stop();
		QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdNonlocalTr call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
		//	// (S. Kazmin) this part do nothing anyway because linkDirs size is zero comment out for future usage
		//	Timer.reset();
		//	Timer.start();
		//	QDPIO::cout << __func__ << ": start AddLinks" << std::endl;
		//	AddLinks(B, F, U, LinkDirs, MaxNLinks, LinkPattern, 0, -1, qio_file, GBB_NLinkPatterns);
		//	Timer.stop();
		//	QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
		//	TotalTime.stop();
		//	QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
}  // end namespace Chroma
