/*! \file npr_vertex_nonlocal_w.cc
 *  \author Stanislav Kazmin (University Leipzig)
 *  \date 2015-11-24
 *  \brief NPR vertex calculations for the non local axial current.
 */
#include <iomanip>      // std::setprecision

#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_nonlocal_w.h"
#include "meas/hadron/npr_vertex_w.h"

namespace Chroma
{

	void BkwdFrwdNonlocal(const LatticePropagator& Prop2, const LatticePropagator& Prop1, const multi1d<LatticeColorMatrix>& U, QDPFileWriter& qio_file, int& GBB_NLinkPatterns, const multi1d<int>& LinkDirs)
	// TODO (S. Kazmin): linkDirs is not used, do we need it? write it because npr routine need it
	{
		// reset time measurement
		StopWatch TotalTime;
		TotalTime.reset();
		TotalTime.start();
		LatticePropagator Prop2_back = Gamma(15) * adj(Prop2) * Gamma(15); // Prop 2 backwards // TODO check adj 
		// calculate one link operator in each direction
		for(int mu = 0; mu < 1; ++mu) // go though all directions // TODO change to all indices with Ns as maximum
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
			// Compute the single site propagator and write it
			DPropagator prop;
			// A_mu (x) = 1/2 * (ψ^adj_x * γ_5 * γ_μ * U_μ * ψ_{x+μ} + ψ^adj_{x+mu} * γ_5 * γ_μ * U^adj_μ * ψ_x)
			// TODO (S. Kazmin): gamma 5 in Prop2 or not? U Adjungation!!!
// 			LatticePropagator tmp = 0.5 * (Prop2 * Gamma(gamma) * U[mu] * shift(Prop1, FORWARD, mu) + shift(Prop2, FORWARD, mu) * Gamma(gamma) * adj(U[mu]) * Prop1);
// 			// The site's worth of data of interest
// 			// sum is over the volume at each site
// 			prop = sum(tmp) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at all sites
			// append to xml file
// 			///////////////// TEST /////////////////////
// 			LatticePropagator tmp1;
// 			LatticePropagator tmp2;
// 			tmp1 = U[mu] * shift(Prop1, FORWARD, mu);
// 			tmp2 = adj(Prop2_back) * Gamma(gamma) * tmp1;
// 			tmp1 = U[mu] * shift(Prop2_back, FORWARD, mu);
// 			tmp2 += adj(tmp1) * Gamma(gamma) * Prop1;
// 			tmp2 = 0.5 * tmp2;
// 			
// 			//tmp = 0.5 * (Prop2 * Gamma(gamma) * (U[mu] * shift(Prop1, FORWARD, mu)) + shift(Prop2, BACKWARD, mu) * Gamma(gamma) * (adj(U[mu]) * Prop1));
// // 			tmp = Prop2 * Gamma(gamma) * (U[mu] * Prop1);
// 			for (int i = 0; i < Layout::vol(); ++i)
// 			{
// 				QDPIO::cout << "\ttrace " << i << ": " << trace(tmp2.elem(i)) << std::endl;
// 			}
// 	// The site's worth of data of interest
// 			prop = sum(tmp2)/Double(Layout::vol()); // and normalize by the volume
// 			QDPIO::cout << "\ttrace mean: " << trace(prop) << std::endl;
// 			///////////////// END TEST /////////////////////
// 			// TEST-OUTPUT //
// 			//int ind = 0;
// 			LatticePropagator Fshift = shift(Prop1, FORWARD, mu);
// 			for (int ind = 0; ind < Layout::vol(); ++ ind)
// 			{
// 			std::cout << std::setprecision(2);
// 			std::cout << "place = " << ind  << "\n";
// 			// F OUTPUT
// 			std::cout << "F = \n";
// 			for(int i=0; i<12; i++)
// 			{
// 				for(int j=0; j<12; j++)
// 				{
// 					std::cout << "(";
// 					std::cout << Prop1.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).real();
// 					std::cout << " ";
// 					std::cout << Prop1.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).imag();
// 					std::cout << ")";
// 				}
// 				std::cout << "\n";
// 			}
// 			// F shift OUTPUT
// 			std::cout << "F shift = \n";
// 			for(int i=0; i<12; i++)
// 			{
// 				for(int j=0; j<12; j++)
// 				{
// 					std::cout << "(";
// 					std::cout << Fshift.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).real();
// 					std::cout << " ";
// 					std::cout << Fshift.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).imag();
// 					std::cout << ")";
// 				}
// 				std::cout << "\n";
// 			}
// 			// B OUTPUT
// 			std::cout << "B = \n";
// 			for(int i=0; i<12; i++)
// 			{
// 				for(int j=0; j<12; j++)
// 				{
// 					std::cout << "(";
// 					std::cout << Prop2_back.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).real();
// 					std::cout << " ";
// 					std::cout << Prop2_back.elem(ind).elem(i % 4, j % 4).elem(i / 4, j / 4).imag();
// 					std::cout << ")";
// 				}
// 				std::cout << "\n";
// 			}
// 			std::cout << "U[ " << mu << "] = \n";
// 			}
// 			//////////////////
			

			pop(record_xml);
			write(qio_file, record_xml, prop);
		}
		// print elapsed time
		TotalTime.stop();
		QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}

//! NPR vertices
	void NprVertexNonlocal(const LatticePropagator& Prop1,const multi1d< LatticeColorMatrix >& U,const unsigned short int MaxNLinks, const BBLinkPattern LinkPattern, QDPFileWriter& qio_file)
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
		// TODO (S. Kazmin): Backward Propagator should be with different source and sink !?
		//LatticePropagator Prop2 = Gamma(15) * adj(Prop1) * Gamma(15); // backward propagator
		BkwdFrwdNonlocal(Prop1, Prop1, U, qio_file, GBB_NLinkPatterns, LinkDirs);
		Timer.stop();
		QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdNonlocalTr call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
//	// (S. Kazmin) this part do nothing anyway because linkDirs size is zero comment out for future usage
//	Timer.reset();
//	Timer.start();
//	QDPIO::cout << __func__ << ": start AddLinks" << std::endl;
//	AddLinks(Prop2, Prop1, U, LinkDirs, MaxNLinks, LinkPattern, 0, -1, qio_file, GBB_NLinkPatterns);
//	Timer.stop();
//	QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
//	TotalTime.stop();
//	QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
}  // end namespace Chroma
