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
		for(int i = 0; i < Ns * Ns; ++i) // go though all gamma operators = 16 possibilities
		{
			XMLBufferWriter record_xml;
			push(record_xml, "Vertex");
			QDPIO::cout << __func__ << ": LinkDirs = " << LinkDirs << "  gamma = " << i << std::endl;
			write(record_xml, "linkDirs", LinkDirs);   // link pattern
			write(record_xml, "gamma", i);
			// counts number of link patterns
			GBB_NLinkPatterns++;
			// Compute the single site propagator and write it
			DPropagator prop;
			// assumes any Gamma5 matrices have already been absorbed into B
			// TODO (S. Kazmin): what does this mean?
			LatticePropagator tmp = B * Gamma(i) * F;
			// The site's worth of data of interest
			// TODO (S. Kazmin): have we to sum? We have only one propagator anyway
			prop = sum(tmp) / Double(Layout::vol()); // and normalize by the volume
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
