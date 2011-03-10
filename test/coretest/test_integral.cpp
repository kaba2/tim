#include "estimation.h"

#include "tim/core/signal_tools.h"
#include "tim/core/delay_embed.h"

#include <pastel/sys/views.h>
#include <pastel/sys/random.h>
#include <pastel/sys/string_algorithms.h>

#include <pastel/math/matrix_tools.h>

#include <pastel/gfx/pcx.h>
#include <pastel/gfx/summedareatable.h>
#include <pastel/gfx/noise.h>

#include <numeric>

using namespace Tim;

namespace
{

	void testIntegralCase(
		integer n, integer k,
		real correlation)
	{
		ENSURE_OP(std::abs(correlation), <, 1);

		const integer m = 2;

		const AlignedBox2 pdfWindow(
			-3, -3, 3, 3);

		const integer xBins = 2048 + 1;
		const integer yBins = xBins;

		const Matrix2 covariance(
			1, correlation,
			correlation, 1);
		const Matrix2 invCovariance = inverse(covariance);

		//const real covDet = determinant(covariance);

		Array<real> pdf(xBins, yBins, 0);

		for (integer y = 0;y < yBins;++y)
		{
			for (integer x = 0;x < xBins;++x)
			{
				const Vector2 uv = 
					pdfWindow.at(Vector2(
					dequantizeUnsigned(x, xBins),
					dequantizeUnsigned(y, yBins)));

				//pdf(x, y) = std::exp(-dot(uv * invCovariance, uv) / 2);
				//pdf(x, y) = perlinNoise(evaluate(uv * 10));
				pdf(x, y) = perlinNoise(evaluate(uv * 3)) * square(perlinNoise(evaluate(uv)));
			}
		}

		const real fNormalization = 
			std::pow((real)k, (real)m) *
			binomial<real>(n - 1, k);

		//log() << "C = " << fNormalization << logNewLine;

		const real pdfSum = std::accumulate(
			pdf.begin(), pdf.end(), (real)0);
		const real binWidth = pdfWindow.extent()[0] / xBins;
		const real binHeight = pdfWindow.extent()[1] / yBins;
		const real binSize = binWidth * binHeight;
		const real invFactor = inverse(pdfSum * binSize);

		for (integer i = 0;i < pdf.size();++i)
		{
			pdf(i) *= invFactor;
		}

		Array<real> summedPdf(pdf.extent());
		
		computeSummedAreaTable(constArrayView(pdf), arrayView(summedPdf));

		ENSURE_OP(
			relativeError<real>(summedPdf(summedPdf.size() - 1) * binSize, 1), <, 0.001);
		
		//Vector2i center(0, yBins - 1);
		Vector2i center(xBins / 2, yBins / 2);

		const integer maxSteps = 
			std::max(
			std::max(center.x() + 1, xBins - center.x()),
			std::max(center.y() + 1, yBins - center.y()));

		Array<real> integ(maxSteps, maxSteps, 0);

		real integral = 0;
		for (integer y = 0;y < maxSteps;++y)
		{
			for (integer x = 0;x < maxSteps;++x)
			{
				const real q00 = 
					summedAreaTable(
					constArrayView(summedPdf), 
					AlignedBox2i(
					center.x() - x, center.y() - y,
					center.x() + x + 1, center.y() + y + 1)) * binSize;

				const real q10 = 
					summedAreaTable(
					constArrayView(summedPdf), 
					AlignedBox2i(
					center.x() - (x + 1), center.y() - y,
					center.x() + (x + 1) + 1, center.y() + y + 1)) * binSize;

				const real q01 = 
					summedAreaTable(
					constArrayView(summedPdf), 
					AlignedBox2i(
					center.x() - x, center.y() - (y + 1),
					center.x() + x + 1, center.y() + (y + 1) + 1)) * binSize;

				const real q11 = 
					summedAreaTable(
					constArrayView(summedPdf), 
					AlignedBox2i(
					center.x() - (x + 1), center.y() - (y + 1),
					center.x() + (x + 1) + 1, center.y() + (y + 1) + 1)) * binSize;

				const real qDiff = 
					std::max(((q11 - q01) - (q10 - q00)) / binSize, (real)0);

				const integer z = std::max(x, y);

				const real qb = 
					std::max(
					summedAreaTable(
					constArrayView(summedPdf), 
					AlignedBox2i(
					center.x() - z, center.y() - z,
					center.x() + z + 1, center.y() + z + 1)) * binSize, (real)0);

				const real q = 
					std::max(q00, (real)0);

				const real f = 
					fNormalization *
					std::pow(q, (real)(k - 1)) * 
					std::pow(1 - qb, (real)(n - k - 1)) * 
					qDiff;

				ENSURE_OP(q, >=, 0);
				ENSURE_OP(qb, >=, 0);

				const real integrand = f * std::log(q);
				integ(x, y) = -integrand;

				/*
				const real integrand = f;
				integ(x, y) = integrand;
				*/

				integral += integrand * binSize;
			}
		}

		//const real proposed = 1;

		const real proposed = digamma<real>(k) - digamma<real>(n) - (real)(m - 1) / k;

		log() << "n = " << n << ", k = " << k << ", cor = " << correlation << logNewLine;
		log() << "  Numeric:  " << integral << logNewLine;
		log() << "  Analytic: " << proposed << logNewLine;
		log() << "  RelError: " << relativeError<real>(integral, proposed) * 100 
			<< "%" << logNewLine;

		const std::string name = 
			"integral_" + integerToString(n) + "_" + integerToString(k, 2) + 
			"_" + realToString(correlation);

		/*
		saveGrayscalePcx(integ, name + "_integrand.pcx", true);
		saveGrayscalePcx(pdf, name + "_pdf.pcx", true);
		saveGrayscalePcx(summedPdf, name + "_cdf.pcx", true);
		*/
	}

	void testIntegral()
	{
		testIntegralCase(10000, 1, 0);
		testIntegralCase(10000, 2, 0);
		testIntegralCase(10000, 4, 0);
		testIntegralCase(10000, 8, 0);
		testIntegralCase(10000, 16, 0);

		testIntegralCase(1000, 1, 0);
		testIntegralCase(2000, 1, 0);
		testIntegralCase(4000, 1, 0);
		testIntegralCase(8000, 1, 0);
		testIntegralCase(16000, 1, 0);

		testIntegralCase(10000, 1, 0.0);
		testIntegralCase(10000, 1, 0.2);
		testIntegralCase(10000, 1, 0.4);
		testIntegralCase(10000, 1, 0.6);
		testIntegralCase(10000, 1, 0.8);
	}

	void addTest()
	{
		timTestList().add("Integral", testIntegral);
	}

	CallFunction run(addTest);

}

