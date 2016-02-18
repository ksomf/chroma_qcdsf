/*! \file npr_vertex_nonlocal_w.cc
 * \author Stanislav Kazmin (University Leipzig)
 * \date 2016-01-15 14:33:00
 * \brief NPR vertex calculations for the non local axial current.
 */

#include "util/ft/sftmom.h"
#include "meas/hadron/npr_vertex_nonlocal_w.h"
#include "meas/hadron/npr_vertex_w.h"

#include <vector>
#include <cmath>

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


		LatticePropagator Diff_A_x = Complex(0.0); 		//// create lattice propagator dA(x)
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
			//// (S. Kazmin): Amu_x_for is forward directed and is placed at x
			//// (S. Kazmin): A_μ(x) = 1/2 *( B(x) γ_μ γ_5 U_μ(x) F(x+μ) + B(x+μ) adj(U_μ(x)) γ_μ γ_5  F(x))
			//// (S. Kazmin): the 1/2 factor is shifted to later calculations to reduce thew number of used factor*matrix operations
			LatticePropagator Amu_x_for =
			adj(B) * Gamma(gamma) * (U[mu] * shift(F, FORWARD, mu))
			+ adj(U[mu] * shift(B, FORWARD, mu)) * Gamma(gamma)  * F;
			//// (S. Kazmin): Amu_x_back is forward directed and is placed at x
			//// (S. Kazmin): A_μ(x-mu) = 1/2 *( B(x-mu) γ_μ γ_5 U_μ(x-mu) F(x) + B(x) adj(U_μ(x-mu)) γ_μ γ_5  F(x-mu))
			//// (S. Kazmin): the 1/2 factor is shifted to later calculations to reduce thew number of used factor*matrix operations
			LatticePropagator Amu_x_back =
			shift(adj(B), BACKWARD, mu) * Gamma(gamma) * (shift(U[mu], BACKWARD, mu) * F)
			+ adj(shift(U[mu], BACKWARD, mu) * B) * Gamma(gamma) * shift(F, BACKWARD, mu);

			//// sum is over the volume at each site
			Amu_x_mean = 0.5 * sum(Amu_x_for) / Double(Layout::vol()); // and normalize by the volume -> mean value of the prop at all sites
			//// append to xml file
			pop(record_xml);
			write(qio_file, record_xml, Amu_x_mean);
			// calculate the dA operator
			Diff_A_x += (Amu_x_for - Amu_x_back);
			// std::cout << "testing " << "\t for  " << Amu_x_for.elem(0).elem(0,0).elem(0,0).real() << std::endl;
			// std::cout << "testing " << "\t back " << Amu_x_back.elem(0).elem(0,0).elem(0,0).real() << std::endl;
			// std::cout << "testing " << "\t diff " << Diff_A_x.elem(0).elem(0,0).elem(0,0).real() << std::endl;
		}
		//// for testing reason add the momentum array manually
		multi1d<int> nrow(4);
		nrow[0] = 4; nrow[1] = 4; nrow[2] = 4; nrow[3] = 8;

		multi1d<int> tmpMom(4);
		std::vector<multi1d<int>> mom;
		// tmpMom[0] = 0; tmpMom[1] = 0; tmpMom[2] = 0; tmpMom[3] = 0;
		// mom.push_back(tmpMom);
		// tmpMom[0] = 0; tmpMom[1] = 0; tmpMom[2] = 0; tmpMom[3] = 1;
		// mom.push_back(tmpMom);
		// tmpMom[0] = 1; tmpMom[1] = 0; tmpMom[2] = 0; tmpMom[3] = 0;
		// mom.push_back(tmpMom);
		// tmpMom[0] = 1; tmpMom[1] = 1; tmpMom[2] = 0; tmpMom[3] = 0;
		// mom.push_back(tmpMom);
		// tmpMom[0] = 1; tmpMom[1] = 1; tmpMom[2] = 1; tmpMom[3] = 0;
		// mom.push_back(tmpMom);
		for (int i = 1; i < 9; ++i)
		{
			tmpMom[0] = i; tmpMom[1] = i; tmpMom[2] = i; tmpMom[3] = 2*i;
			mom.push_back(tmpMom);
		}
		// for (auto iter : mom)
		// {
		// 	std::cout << iter[0] << iter[1] << iter[2] << iter[3] << "\n";
		// }
		//// calculate the sin^2(p) values
		std::vector<double> sin2p;
		for (auto iter : mom)
		{
			double sin2pValue = 0.0;
			for (int mu = 0; mu < Nd; ++mu)
			{
				sin2pValue += pow(sin(iter[mu] * 2. * M_PI / (double)(nrow[mu])), 2.0);
			}
			sin2p.push_back(sin2pValue);
			// std::cout << iter[0] << iter[1] << iter[2] << iter[3] << " : sin^2(p) = ";
			// std::cout << sin2pValue << "\n";
		}
		//// go through all momenta
		for (auto iter : mom)
		{
			
		}




		// // check the matrix
		// for(int i = 0; i < Layout::vol(); ++i)
		// {
		// 	std::cout << "site " << i << "\t Diff_A_x(0,0) = " << Diff_A_x.elem(i).elem(0,0).elem(0,0).real() << std::endl;
		// }
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
