/*! \file npr_vertex_nonlocal_w.cc
 * \author Stanislav Kazmin (University Leipzig)
 * \date 2016-01-15 14:33:00
 * \brief NPR vertex calculations for the non local axial current.
 */

#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_nonlocal_w.h"
#include "meas/hadron/npr_vertex_w.h"

namespace Chroma
{
	
	void BkwdFrwdNonlocal(const LatticePropagator& F, const multi1d<LatticeColorMatrix>& U, QDPFileWriter& qio_file, int& GBB_NLinkPatterns, const multi1d<int>& LinkDirs)
	//// (S. Kazmin): linkDirs is not used, do we need it? write it because npr routine need it
	{
		// reset time measurement
		StopWatch TotalTime;
		TotalTime.reset();
		TotalTime.start();
		LatticePropagator B = Gamma(15) * F * Gamma(15); //// (S. Kazmin): backwards propagator NOT adj.
		// calculate one link operator in each direction
		for(int mu = 0; mu < Nd; ++mu) // go though all directions
		{
			int gamma = 1 << mu; // get gamma value from bits-shift of 1.
			gamma = 15 - gamma; // get the gamma index of the according γ_μ γ_5 operator
			XMLBufferWriter record_xml;
			push(record_xml, "Vertex");
			QDPIO::cout << __func__ << " calculate point split operator A_mu  for gamma = " << gamma << std::endl;
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
			//// The site's worth of data of interest
			//// sum is over the volume at each site
			//// TODO (S. Kazmin): will try to implement the full local divergence and ZA calculation here for test reasons
			Amu_x_mean = 0.5 * sum(Amu_x) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at all sites
			//// append to xml file
			pop(record_xml);
			write(qio_file, record_xml, Amu_x_mean);
		}
		// print elapsed time
		TotalTime.stop();
		QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
	//! NPR vertices
	void NprVertexNonlocal(const LatticePropagator & F, const multi1d< LatticeColorMatrix >& U, const unsigned short int MaxNLinks, const BBLinkPattern LinkPattern, QDPFileWriter & qio_file)
	{
		StopWatch Timer;
		Timer.reset();
		Timer.start();
		int GBB_NLinkPatterns;
		// calculate building blocks
		QDPIO::cout << __func__ << ": start BkwdFrwdNonlocal (split point operator A_mu)" << std::endl;
		const int NLinks = 0;
		multi1d<int> LinkDirs(0);
		//// (S. Kazmin): Backward propagator is calculated in BkwdFrwdNonlocal function
		BkwdFrwdNonlocal(F, U, qio_file, GBB_NLinkPatterns, LinkDirs);
		Timer.stop();
		QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdNonlocal call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
		return;
	}
}  // end namespace Chroma
