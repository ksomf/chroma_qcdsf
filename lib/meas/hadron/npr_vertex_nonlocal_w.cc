/*! \file npr_vertex_nonlocal_w.cc
 *  \author Stanislav Kazmin (University Leipzig)
 *  \date 2015-11-24
 *  \brief NPR vertex calculations for the non local axial current.
 */

#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_nonlocal_w.h"
#include "meas/hadron/npr_vertex_w.h"

namespace Chroma
{

	void BkwdFrwdNonlocal(const LatticePropagator& B,
						  const LatticePropagator& F,
						  const multi1d<LatticeColorMatrix>& U,
						  QDPFileWriter& qio_file,
						  int& GBB_NLinkPatterns,
						  const multi1d<int>& LinkDirs)// TODO (S. Kazmin): linkDirs is not used, do we need it? write it because npr routine need it
	{
		// reset time measurement
		StopWatch TotalTime;
		TotalTime.reset();
		TotalTime.start();
		// calculate one link operator
		for(int mu = 0; mu < Ns; ++mu) // go though all directions
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
			// TODO (S. Kazmin): gamma 5 in B or not? U Adjungation!!!
			LatticePropagator tmp = 0.5 * (B * Gamma(gamma) * U[mu] * shift(F, FORWARD, mu) + shift(B, FORWARD, mu) * Gamma(gamma) * adj(U[mu]) * F);
			// The site's worth of data of interest
			// sum is over the volume at each site
			prop = sum(tmp) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at at sites
			// append to xml file
			pop(record_xml);
			write(qio_file, record_xml, prop);
		}
		// print elapsed time
		TotalTime.stop();
		QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
//! NPR vertices
	void NprVertexNonlocal(const LatticePropagator& F,
						   const multi1d< LatticeColorMatrix >& U,
						   const unsigned short int MaxNLinks,
						   const BBLinkPattern LinkPattern,
						   QDPFileWriter& qio_file)
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
		// TODO (S. Kazmin): gamm15 or gamma0?
		LatticePropagator B = Gamma(15) * adj(F) * Gamma(15); // backward propagator
		BkwdFrwdNonlocal(B, F, U, qio_file, GBB_NLinkPatterns, LinkDirs);
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
