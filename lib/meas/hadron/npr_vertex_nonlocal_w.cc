/*! \file npr_vertex_nonlocal_w.cc
 * \author Stanislav Kazmin (University Leipzig)
 * \date 2016-01-13 10:44:30
 *  \brief NPR vertex calculations for the non local axial current.
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
			DPropagator Amu_x_Mean;
			//// (S. Kazmin): Amu_x is forward directed and is placed at x
			//// (S. Kazmin): A_μ(x) = 1/2 *( B γ_μ γ_5 U_μ F_μ + B_μ adj(U_μ) γ_μ γ_5  F)
			LatticePropagator Amu_x = 0.5 * (adj(B) * Gamma(gamma) * (U[mu] * shift(F, FORWARD, mu)) + adj(U[mu] * shift(B, FORWARD, mu)) * Gamma(gamma)  * F);
			LatticePropagator Amu_xMinusMu = 0.5 * (adj(B) * Gamma(gamma) * (U[mu] * shift(F, FORWARD, mu)) + adj(U[mu] * shift(B, FORWARD, mu)) * Gamma(gamma)  * F);
			// The site's worth of data of interest
			// sum is over the volume at each site
			Amu_x_Mean = sum(Amu_x) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at all sites
			// append to xml file
			pop(record_xml);
			write(qio_file, record_xml, Amu_x_Mean);
		}
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
