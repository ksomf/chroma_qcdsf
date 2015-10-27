// $Id: slrc_feynhell_fermact_params_w.cc,v 3.4 2008-11-10 17:59:07 bjoo Exp $
/*! \file
 *  \brief Clover fermion action parameters
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/slrc_feynhell_fermact_params_w.h"

#include "io/param_io.h"
#include "util/ft/sftmom.h"

namespace Chroma
{

	SLRCFeynHellFermActParams::SLRCFeynHellFermActParams()
	{
		QDPIO::cout << "Creating SLRCFeynHellFermActParams" << std::endl;
	}

	SLRCFeynHellFermActParams::SLRCFeynHellFermActParams(XMLReader& xml, const std::string& path)
	{
		// XMLReader paramtop(xml, path);

		// numParam = paramtop.count("FeynHellParam/elem");
		// QDPIO::cout << "In constructor, xml path is: " << path << std::endl;
		// QDPIO::cout << "Number of FeynHellParam elems is: "
		// 	    << numParam << std::endl << std::flush;
		// if(numParam == 0)
		// {
		// 	QDPIO::cerr << "ERROR: No FeynHellParam elems found"
		// 		    << std::endl << std::flush;
		// 	QDP_abort(1);
		// }
		// read(paramtop, "FeynHellParam", FHparam);
	  read(xml, path, *this);

	}

	void read(XMLReader& xml, const std::string& path, SLRCFeynHellFermActParams& param)
	{
		XMLReader paramtop(xml,path);
		param.numParam = paramtop.count("FeynHellParam/elem");
		QDPIO::cout << "In read, xml path is: " << path << std::endl;
		QDPIO::cout << "Number of FeynHellParam elems is: "
			    << param.numParam << std::endl << std::flush;
		if(param.numParam == 0)
		{
			QDPIO::cerr << "ERROR: No FeynHellParam elems found " << std::endl << std::flush;
			QDP_abort(1);
		}
		read(paramtop, "FeynHellParam", param.FHparam);
	}

	void read(XMLReader& xml, const std::string& path, SLRCFeynHellFermActParams::FHParam& param)
	{
		XMLReader paramtop(xml, path);

		Real lambdaReal;
		Real lambdaImag;

		if (paramtop.count("LambdaReal") != 0)
			read(paramtop, "LambdaReal", lambdaReal);
		else
		{
			QDPIO::cerr << "Real part of lambda not given" << std::endl << std::flush;
			QDP_abort(1);
		}
		if (paramtop.count("LambdaImag") != 0)
			read(paramtop, "LambdaImag", lambdaImag);
		else
		{
			QDPIO::cerr << "Imaginary part of lambda not given" << std::endl << std::flush;
			QDP_abort(1);
		}
		param.lambda = cmplx(lambdaReal,lambdaImag);

		Real noise_real;
		Real noise_imag;

		if (paramtop.count("NoiseReal") != 0 && paramtop.count("NoiseImag") != 0)
		{
			read(paramtop, "NoiseReal", noise_real);
			read(paramtop, "NoiseImag", noise_imag);
			QDPIO::cout << "Found noise" << std::endl << std::flush;
		}
		else
		{
			noise_real = 1.0;
			noise_imag = 0.0;
			QDPIO::cout << "No noise found (or only one of real/imag part given)" << std::endl << std::flush;
		}
		param.noise = cmplx(noise_real,noise_imag);

		if (paramtop.count("Operator") != 0)
			read(paramtop, "Operator", param.op);
		else
		{
			QDPIO::cerr << "Operator not given" << std::endl << std::flush;
			QDP_abort(1);
		}

		if (paramtop.count("Source") != 0)
			read(paramtop, "Source", param.source);
		else
		{
			QDPIO::cerr << "Source not given" << std::endl << std::flush;
			QDP_abort(1);
		}

		if (paramtop.count("Momentum") != 0)
			read(paramtop, "Momentum", param.mom);
		else
		{
			QDPIO::cerr << "Momentum not given" << std::endl << std::flush;
			QDP_abort(1);
		}
		param.phases = SftMom(0,param.source,param.mom,false,Nd-1)[0];

		//Real noiseReal;
		//Real noiseImag;
		//read(paramtop, "NoiseReal", noiseReal);
		//read(paramtop, "NoiseImag", noiseImag);
		//param.noise = cmplx(noiseReal,noiseImag);
		//
		QDPIO::cout << "Lambda is " << param.lambda << " with op: " << param.op << std::endl << std::flush;
		QDPIO::cout << "Momentum is " << param.mom[0] << ", " << param.mom[1] << ", " << param.mom[2] << std::endl << std::flush;
		QDPIO::cout << "Noise is " << param.noise << std::endl << std::flush;

	}

	void write(XMLWriter& xml, const std::string& path, const SLRCFeynHellFermActParams& param)
	{
		push(xml, path);

		write(xml, "FeynHellParam", param.FHparam);

		pop(xml);
	}
	void write(XMLWriter& xml, const std::string& path, const SLRCFeynHellFermActParams::FHParam& param)
	{
		push(xml, path);

		write(xml, "Lambda", param.lambda);
		write(xml, "Operator", param.op);
		write(xml, "Momentum", param.mom);
		write(xml, "Noise", param.noise);
		write(xml, "Source", param.source);

		pop(xml);
	}


}
