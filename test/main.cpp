

#include "zeta.h"
#include <vector>
#include <iostream>

std::vector<std::vector<numericalzeta::arg_t>> resultMatrix;

#ifdef _MSC_VER
	#if _MSC_VER >= 1914
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	constexpr
	#else
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	const					// MSVC v19.10 seems to not support constexpr complex<long double>::complex(long double) (MSVC v19.14 & up seems to work fine, as well as CLANG 6.0.0 & up and GCC 6.1 & up)
	#endif
#else
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	constexpr
#endif


_NOT_CONSTEXPR_IFF_UNSUPPORTED numericalzeta::arg_t zeta_1p1i		= { 0.58215805975200364819946316791425920187798931682653464572,		0.92684856433080707653642431391750077405345489387394342768 };
_NOT_CONSTEXPR_IFF_UNSUPPORTED numericalzeta::arg_t zeta_100p100i	= { 1.00000000000000000000000000000077318636991549092654,			-1.56474806799752257276050601547768116560141375941075e-31 };
_NOT_CONSTEXPR_IFF_UNSUPPORTED numericalzeta::arg_t zeta_05p4000i	= { 0.0172858,														-0.0446159 };
_NOT_CONSTEXPR_IFF_UNSUPPORTED numericalzeta::arg_t zeta_05p100i	= { 2.69262,														-0.0203860 };

int main()
{
	numericalzeta::arg_t z, zeta_z, zeta_z_def, zeta_z_cvz;
	uint32_t precision, iterations, itrs_needed_def, itrs_needed_cvz;

	z = { 100.0, 0 };			zeta_z = { 0 };					precision = 10;					iterations = 100;
	itrs_needed_def = numericalzeta::RiemannZetaFnDEFitrsNeeded(z, precision);		itrs_needed_cvz = numericalzeta::RiemannZetaFnCVZitrsNeeded(z, precision);
	zeta_z_def = numericalzeta::CalculateRiemannZetaFnDEF(z, iterations);					zeta_z_cvz = numericalzeta::CalculateRiemannZetaFnCVZ(z, iterations);

	std::cout << "[" << precision << "/" << iterations << "]\t" << z << "\tbyDef: [" << itrs_needed_def << "]\t" << zeta_z_def << "\tbyCVZ: [" << itrs_needed_cvz << "]\t" << zeta_z_cvz << "\tTruth:" << zeta_z << "\n";

	// RiemannZetaFnDEFitrsNeeded not working from here!

	z = { 100.0,100.0 };			zeta_z = zeta_100p100i;				precision = 10;					iterations = 10;
	itrs_needed_def = numericalzeta::RiemannZetaFnDEFitrsNeeded(z, precision);		itrs_needed_cvz = numericalzeta::RiemannZetaFnCVZitrsNeeded(z, precision);
	zeta_z_def = numericalzeta::CalculateRiemannZetaFnDEF(z, iterations);			zeta_z_cvz = numericalzeta::CalculateRiemannZetaFnCVZ(z, iterations);

	std::cout << "[" << precision << "/" << iterations << "]\t" << z << "\tbyDef: [" << itrs_needed_def << "]\t" << zeta_z_def << "\tbyCVZ: [" << itrs_needed_cvz << "]\t" << zeta_z_cvz << "\tTruth:" << zeta_z << "\n";

	z = { 1.0,1.0 };				zeta_z = zeta_1p1i;					precision = 10;					iterations = 10;
	itrs_needed_def = numericalzeta::RiemannZetaFnDEFitrsNeeded(z, precision);		itrs_needed_cvz = numericalzeta::RiemannZetaFnCVZitrsNeeded(z, precision);
	zeta_z_def = numericalzeta::CalculateRiemannZetaFnDEF(z, iterations);			zeta_z_cvz = numericalzeta::CalculateRiemannZetaFnCVZ(z, iterations);

	std::cout << "[" << precision << "/" << iterations << "]\t" << z << "\tbyDef: [" << itrs_needed_def << "]\t" << zeta_z_def << "\tbyCVZ: [" << itrs_needed_cvz << "]\t" << zeta_z_cvz << "\tTruth:" << zeta_z << "\n";

	z = { 0.5,100.0 };			zeta_z = zeta_05p100i;				precision = 10;					iterations = 100;
	itrs_needed_def = numericalzeta::RiemannZetaFnDEFitrsNeeded(z, precision);		itrs_needed_cvz = numericalzeta::RiemannZetaFnCVZitrsNeeded(z, precision);
	zeta_z_def = numericalzeta::CalculateRiemannZetaFnDEF(z, 100000);			zeta_z_cvz = numericalzeta::CalculateRiemannZetaFnCVZ(z, iterations);

	std::cout << "[" << precision << "/" << iterations << "]\t" << z << "\tbyDef: [" << itrs_needed_def << "]\t" << zeta_z_def << "\tbyCVZ: [" << itrs_needed_cvz << "]\t" << zeta_z_cvz << "\tTruth:" << zeta_z << "\n";

	z = { 100.0, 0.5 };			zeta_z = { 0 };					precision = 10;					iterations = 100;
	itrs_needed_def = numericalzeta::RiemannZetaFnDEFitrsNeeded(z, precision);		itrs_needed_cvz = numericalzeta::RiemannZetaFnCVZitrsNeeded(z, precision);
	zeta_z_def = numericalzeta::CalculateRiemannZetaFnDEF(z, iterations);					zeta_z_cvz = numericalzeta::CalculateRiemannZetaFnCVZ(z, iterations);

	std::cout << "[" << precision << "/" << iterations << "]\t" << z << "\tbyDef: [" << itrs_needed_def << "]\t" << zeta_z_def << "\tbyCVZ: [" << itrs_needed_cvz << "]\t" << zeta_z_cvz << "\tTruth:" << zeta_z << "\n";

	// testing main function

	numericalzeta::value_ready_func_t populatorFunction = numericalzeta::ValueReadyFuncFactoryWithResultMatrix({ 0.3,0.0 }, { 0.6,100.0 }, resultMatrix, { 0.1, 0.1 });
	numericalzeta::CalculateRiemannZetaFn({ 0.3,0.0 }, { 0.6,100.0 }, populatorFunction, { 0.1, 0.1 });

	std::cout << "resultMatrix: " << resultMatrix.size() << "x" << (resultMatrix.empty() ? 0 : resultMatrix[0].size()) << "\n";

	for (const auto& row : resultMatrix)
	{
		for (const auto& elem : row)
		{
			std::cout << elem << "\t";
		}
		std::cout << std::endl;
	}

}
