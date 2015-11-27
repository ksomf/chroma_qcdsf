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
					  QDPFileWriter& qio_file,
					  int& GBB_NLinkPatterns,
					  const multi1d< int >& LinkDirs)
{
	StopWatch TotalTime;
	TotalTime.reset();
	TotalTime.start();
	for(int i = 0; i < Ns * Ns; i ++)
	{
		XMLBufferWriter record_xml;
		push(record_xml, "Vertex");
		QDPIO::cout << __func__ << ": LinkDirs = " << LinkDirs << "  gamma = " << i << std::endl;
		write(record_xml, "linkDirs", LinkDirs);   // link pattern
		write(record_xml, "gamma", i);
		// counts number of link patterns
		GBB_NLinkPatterns++;
		// assumes any Gamma5 matrices have already been absorbed
		int G5 = Ns * Ns - 1;
		// Compute the single site propagator and write it
		DPropagator prop;
		{
			// assumes any Gamma5 matrices have already been absorbed into B
			LatticePropagator tmp = B * Gamma(i) * F;
			// The site's worth of data of interest
			prop = sum(tmp) / Double(Layout::vol()); // and normalize by the volume
		}
		pop(record_xml);
		write(qio_file, record_xml, prop);
	}
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
	StopWatch TotalTime;
	TotalTime.reset();
	TotalTime.start();
	StopWatch Timer;
	int GBB_NLinkPatterns;
	//#################################################################################//
	// open building blocks data files                                                 //
	//#################################################################################//
	Timer.reset();
	Timer.start();
	//#################################################################################//
	// calculate building blocks                                                       //
	//#################################################################################//
	QDPIO::cout << __func__ << ": start BkwdFrwdNonlocal" << std::endl;
	const int NLinks = 0;
	multi1d< int > LinkDirs(0);
	LatticePropagator B = Gamma(15) * adj(F) * Gamma(15);
	BkwdFrwdNonlocal(B, F, qio_file, GBB_NLinkPatterns, LinkDirs);
	Timer.stop();
	QDPIO::cout << __func__ << ": total time for 0 links (single BkwdFrwdNonlocalTr call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
	Timer.reset();
	Timer.start();
	QDPIO::cout << __func__ << ": start AddLinks" << std::endl;
	AddLinks(B, F, U, LinkDirs, MaxNLinks, LinkPattern, 0, -1, qio_file, GBB_NLinkPatterns);
	Timer.stop();
	QDPIO::cout << __func__ << ": total time for remaining links (outermost AddLinks call) = " << Timer.getTimeInSeconds() << " seconds" << std::endl;
	TotalTime.stop();
	QDPIO::cout << __func__ << ": total time = " << TotalTime.getTimeInSeconds() << " seconds" << std::endl;
	return;
}

}  // end namespace Chroma
