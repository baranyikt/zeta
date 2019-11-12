
#define ZETA_DEFINITION_USE_REVERSE_ITERATION

#ifdef _MSC_VER
	#if _MSC_VER >= 1914
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	constexpr
	#else
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	const					// MSVC v19.10 seems to not support constexpr complex<long double>::complex(long double) (MSVC v19.14 & up seems to work fine, as well as CLANG 6.0.0 & up and GCC 6.1 & up)
	#endif
#else
		#define _NOT_CONSTEXPR_IFF_UNSUPPORTED	constexpr
#endif

#include "zeta.h"

#define __MATH_LONG_DOUBLE_CONSTANTS
#include <cmath>

// check for availability of long double PI constant (sometimes defined in cmath header)
#ifndef M_PIl
// we need predefined PI constant 
#define	M_PIl	3.141592653589793238462643383279502884L
#endif

namespace numericalzeta {

inline void PrepareResultMatrix(const arg_t& leftEndPoint, const arg_t& rightEndPoint, std::vector<std::vector<arg_t>>& mtx, size_t& realdim, size_t& imagdim, arg_t step = DefaultStep)
{
	realdim = ceil((real(rightEndPoint) - real(leftEndPoint)) / real(step)) + 1;
	imagdim = ceil((imag(rightEndPoint) - imag(leftEndPoint)) / imag(step)) + 1;
	mtx.resize(realdim);
	for (auto& vec : mtx)
		vec.resize(imagdim);
}

value_ready_func_t ValueReadyFuncFactoryWithResultMatrix(arg_t leftEndPoint, arg_t rightEndPoint, std::vector<std::vector<arg_t>>& mtx, arg_t step)
{
	size_t realdim, imagdim;
	PrepareResultMatrix(leftEndPoint, rightEndPoint, mtx, realdim, imagdim, step);
	auto retval = [leftEndPoint, rightEndPoint, realdim, imagdim, step, &mtx](const arg_t& z, const arg_t& zeta_z) {
		const size_t realpos = 
			fabs(real(z) - real(leftEndPoint)) < FLOATINGPOINT_EQ_TOLERANCE ? 
			0 :
			((real(z) - real(leftEndPoint)) / real(step)) + 0.5;
		const size_t imagpos = 
			fabs(imag(z) - imag(leftEndPoint) < FLOATINGPOINT_EQ_TOLERANCE ?
			0 : 
			((imag(z) - imag(leftEndPoint)) / imag(step))) + 0.5;
		mtx[realpos][imagpos] = zeta_z;
	};
	return retval;
}

intrmed_int_t RiemannZetaFnCVZitrsNeeded(arg_t z, intrmed_int_t decimal_precision)
{
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _1_3 = static_cast<arg_base_t>(1.3L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _0_9 = static_cast<arg_base_t>(0.9L);
	auto _dec_prec_casted = static_cast<arg_base_t>(decimal_precision);
	return static_cast<intrmed_int_t>(ceil(_1_3 * _dec_prec_casted + _0_9*abs(z.imag())));
}

inline intrmed_t CVZseries(intrmed_int_t k, intrmed_t z)
{
	_NOT_CONSTEXPR_IFF_UNSUPPORTED intrmed_t _1_0(1.0L);
	const auto _complex_k = static_cast<intrmed_t>(k);;
	return _1_0 / pow(_complex_k + _1_0, z);
}

arg_t CalculateRiemannZetaFnCVZ(arg_t z, intrmed_int_t iterations)
{
	// Cohen-Villegas-Zagier method (for smaller z values)
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _0_5 = static_cast<intrmed_t>(0.5L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _1_0 = static_cast<intrmed_t>(1.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _2_0 = static_cast<intrmed_t>(2.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _3_0 = static_cast<intrmed_t>(3.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _8_0 = static_cast<intrmed_t>(8.0L);
	const auto _z = static_cast<intrmed_t>(z);

	// CVZ step #1 -- calculate Eta(z) = Sum_0^infinity (-1)^n / (n+1)^z using Convergence Acceleration of Alternating Series paper
	const intrmed_t d_proto = pow(_3_0 + sqrt(_8_0), static_cast<intrmed_t>(iterations));
	const intrmed_t d_init = (d_proto + _1_0 / d_proto) * _0_5;
	intrmed_t d = d_init;
	intrmed_t c = -d_init;
	intrmed_t b = -1;
	intrmed_t s = 0;
	for (intrmed_int_t k = 0; k < iterations; ++k)
	{
		const intrmed_t a_k = CVZseries(k, _z);
		c = b - c;
		s = s + c * a_k;

		const auto _k		= static_cast<intrmed_t>(k);
		const auto _itrs	= static_cast<intrmed_t>(iterations);

		b *= (_k + _itrs) * (_k - _itrs);
		b /= (_k + _0_5) * (_k + _1_0);
	}
	// CVZ step #2 -- use (1 - 2^(1-z)) Zeta(z) = Eta(z)
	const intrmed_t eta_z = s / d;
	const intrmed_t zeta_z = eta_z / (_1_0 - pow(_2_0, _1_0 - _z));
	return static_cast<arg_t>(zeta_z);
}

/// prerequisites: real(z) > 0, imag(z) == 0!
intrmed_int_t RiemannZetaFnDEFitrsNeeded(arg_t z, intrmed_int_t decimal_precision)
{
	if (real(z) <= 0 || imag(z) != 0)
		return std::numeric_limits<intrmed_int_t>::max();
	// or it could be:
	// if (real(z) <= 0 || imag(z) != 0) throw std::invalid_argument("RiemannZetaFnDEFitrsNeeded: value out of range");
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _1_0	= static_cast<arg_base_t>(1.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _2_0	= static_cast<arg_base_t>(2.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _10_0	= static_cast<arg_base_t>(10.0L);
	auto _dec_prec_casted = static_cast<arg_base_t>(decimal_precision);

	return ceil(pow(_10_0,  (-_dec_prec_casted*(_1_0-real(z))) / (pow(_1_0-z.real(),_2_0) + pow(z.imag(),_2_0))  ));
}

arg_t CalculateRiemannZetaFnDEF(arg_t z, intrmed_int_t iterations)
{
	auto retval = static_cast<intrmed_t>( 0.0L );
	const auto z_casted = static_cast<intrmed_t>(z);
#ifdef ZETA_DEFINITION_USE_REVERSE_ITERATION
	for (intrmed_int_t j = iterations+1; j >= 1; --j)
#else
	for (intrmed_int_t j = 1; j < iterations + 1; ++j)
#endif
	{
		const auto j_casted = static_cast<intrmed_t>(j);
		retval += pow(j_casted,-z_casted);
	}
	return static_cast<arg_t>(retval);
}

inline arg_t drop_imag(intrmed_t z) 
{ 
	if (fabs(z.imag()) < FLOATINGPOINT_EQ_TOLERANCE) 
		return{ static_cast<arg_base_t>(z.real()), 0 }; 
	else 
		return static_cast<arg_t>(z); 
}

/// Fixed-precision gamma from https://en.wikipedia.org/wiki/Lanczos_approximation
inline arg_t gamma(arg_t z)
{
	_NOT_CONSTEXPR_IFF_UNSUPPORTED intrmed_t p[8] = {	676.5203681218851,
														-1259.1392167224028,
														771.32342877765313,
														-176.61502916214059,
														12.507343278686905,
														-0.13857109526572012,
														9.9843695780195716e-6,
														1.5056327351493116e-7		};
	_NOT_CONSTEXPR_IFF_UNSUPPORTED arg_t xs = 			0.99999999999980993;
	constexpr size_t _psize = sizeof(p) / sizeof(p[0]);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _0_5 	= static_cast<intrmed_t>(0.5L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _1_0 	= static_cast<intrmed_t>(1.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _pi 	= static_cast<intrmed_t>(M_PIl);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _psize_casted 	= static_cast<intrmed_t>(_psize);
	const auto _sqrt_two_pi = static_cast<intrmed_t>(sqrt(2.0L*_pi));

	auto z_casted = static_cast<intrmed_t>(z);
	if (z.real() < static_cast<arg_base_t>(0.5L))
	{
		return static_cast<arg_t>(drop_imag(_pi / (sin(_pi * z_casted) * static_cast<intrmed_t>(gamma(static_cast<arg_t>(_1_0 - z_casted))))));
	}
	else
	{
        z_casted -= 1;
        intrmed_t x = xs;
		for (size_t idx = 0; idx < _psize; ++idx)
		{
			auto idx_casted = static_cast<intrmed_t>(idx);
			x += p[idx] / (z_casted + idx_casted + _1_0);
		}
		intrmed_t t = z_casted + _psize_casted - _0_5;
		intrmed_t y = _sqrt_two_pi * pow(t,z_casted + _0_5) * exp(-t) * x;
		return drop_imag(y);
	}
}

inline void CalculateRiemannZetaFnReflectionFormula(arg_t z, arg_t zeta_z, arg_t& reflected_z, arg_t& reflected_zeta_z)
{
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _1_0		= static_cast<arg_base_t>(1.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _2_0		= static_cast<arg_base_t>(2.0L);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _two_pi		= static_cast<arg_base_t>(2.0L * M_PIl);
	_NOT_CONSTEXPR_IFF_UNSUPPORTED auto _pi_half	= static_cast<arg_base_t>(0.5L * M_PIl);
	reflected_z			= _1_0 - z;
	reflected_zeta_z	= _2_0/pow(_two_pi, z) * cos(z * _pi_half) * gamma(z) * zeta_z;
}

inline void CalculateRiemannZetaFnConjugateFormula(arg_t z, arg_t zeta_z, arg_t& conj_z, arg_t& conj_zeta_z)
{
	conj_z		= conj(z);
	conj_zeta_z	= conj(zeta_z);
}

/*
Im z < 0
	conjugate formula
	zeta(conj(z)) = conj(zeta(z))
Re z <= 0
	reflection formula
	Zeta(1-z) = 2/(2 pi)^z * cos (z pi / 2) Gamma(z) Zeta(z)
				|----------------------------------|
								1/chi(z)
*/

void CalculateRiemannZetaFn(
	arg_t leftEndPoint, 
	arg_t rightEndPoint, 
	value_ready_func_t valuepair_ready, 
	arg_t step, 
	uint32_t decimal_precision, 
	progress_func_t progress)
{
	using std::swap;
	if (real(leftEndPoint) > real(rightEndPoint) || imag(leftEndPoint) > imag(rightEndPoint))
		return;
	
	_NOT_CONSTEXPR_IFF_UNSUPPORTED arg_t reflection_center{ 1.0L, 0.0L };					// 1 + 0i: the point to which we need to reflect the requested region 
																	// for the conjugate & reflection formula to work (e.g. to calculate 
																	// zeta(-1+0i), we'll need the value of zeta(1-(-1+0i)) = zeta(2) for 
																	// the reflection formula)
	_NOT_CONSTEXPR_IFF_UNSUPPORTED arg_t direct_methods_working_endpoint { 0.0, 0.0 };	// 0 + 0i: right and upwards from this, we've got direct methods to
																	// calculate zeta(z), while left or down from this we'll need 
																	// reflection and conjugate formulas

	// all these bools indicate special cases that are speeding up the last step of calculation
	const bool direct_methods_enough = real(leftEndPoint) >= real(direct_methods_working_endpoint) && imag(leftEndPoint) >= imag(direct_methods_working_endpoint);
	const bool symmetric_to_center_real = ((real(leftEndPoint) + real(rightEndPoint)) / 2.0 - real(reflection_center)) < FLOATINGPOINT_EQ_TOLERANCE;
	const bool symmetric_to_center_imag = ((imag(leftEndPoint) + imag(rightEndPoint)) / 2.0 - imag(reflection_center)) < FLOATINGPOINT_EQ_TOLERANCE;
	const bool symmetric_to_center_both = symmetric_to_center_real && symmetric_to_center_imag;

	// -- Reflection and extension of original region --
	//
	// We will use direct methods for region { Re z >= 0, Im z >= 0 }. Outside that, we'll use reflection formula for values { Re z < 0 },
	// but those will need the results of the direct methods applied for the reflection of the original region to the { Re z = 1 } line.
	// Thus we'll need to reflect the requested region to that line, but only the negative part of it. (In other words, on the 
	// { 0 <= Im Z <= 1 } critical stripe, we could calculate zeta values in both ways: not only by using direct methods on the original 
	// region, but also by applying those methods to the reflected { 1 <= Im z <= 2 } region. We'll stick with the first one in this case.)
	// The case is similar for { Im z < 0 } values, only a bit simpler, because the reflection line and the boundary line of scope
	// of direct methods are the same { Im z = 0 }. (There is no such two-way calculable region as in real reflection case.) Below
	// that line, we'll need conjugate formula and to use that, we'll need zeta values for the reflection of the original region to
	// the same { Im z = 0 } line.
	// Thus, we divide the complex plane by the real line and by line { Re z = 1 }, we reflect the originally requested region to both lines,
	// take the union of the resulted sets, and keep the part that is approachable by direct methods in the end. This way all the direct 
	// calculations needed either originally or by the reflection/conjugate formulas will be done.
	// [roiLeft,roiRight] rectangle will be the extended region of interest intersection with 
	// { a + ib | a >= 0, b >= 0 }, so using notation anything^w for anything reflected to point w:
	//		[roiLeft,roiRight] = 	[leftEndPoint,rightEndPoint] 
	//									union 
	//								[rightEndPoint^w,leftEndPoint^w]
	//									intersect
	//								{ a + ib | a >= 0, b >= 0 }
	// (This whole process is of course unncessary if the originally requested region is entirely in { a + ib | a >= 0, b >= 0 }.)
	//
	// Reflections can be done dimensionwise.
	arg_t roiLeft, roiRight;
	if (!direct_methods_enough)
	{
		// Reflecting both endpoints horizontally to the right side of reflection_center point
		arg_base_t real_leftEndPoint_xferred_to_quarter1 = real(leftEndPoint)		>= real(reflection_center) ? real(leftEndPoint)		: 2*real(reflection_center)-real(leftEndPoint);
		arg_base_t real_rightEndPoint_xferred_to_quarter1 = real(rightEndPoint)		>= real(reflection_center) ? real(rightEndPoint)	: 2*real(reflection_center)-real(rightEndPoint);

		// and vertically to the upper side of reflection_center point
		arg_base_t imag_leftEndPoint_xferred_to_quarter1	= imag(leftEndPoint)	>= imag(reflection_center) ? imag(leftEndPoint)	: 2*imag(reflection_center)-imag(leftEndPoint);
		arg_base_t imag_rightEndPoint_xferred_to_quarter1	= imag(rightEndPoint)	>= imag(reflection_center) ? imag(rightEndPoint)	: 2*imag(reflection_center)-imag(rightEndPoint);

		// This could have changed their sort order, switch it back if needed
		if (real_leftEndPoint_xferred_to_quarter1 > real_rightEndPoint_xferred_to_quarter1)
			swap(real_leftEndPoint_xferred_to_quarter1, real_rightEndPoint_xferred_to_quarter1);
		if (imag_leftEndPoint_xferred_to_quarter1 > imag_rightEndPoint_xferred_to_quarter1)
			swap(imag_leftEndPoint_xferred_to_quarter1, imag_rightEndPoint_xferred_to_quarter1);
		roiLeft = { real_leftEndPoint_xferred_to_quarter1, imag_leftEndPoint_xferred_to_quarter1 };
		roiRight = { real_rightEndPoint_xferred_to_quarter1, imag_rightEndPoint_xferred_to_quarter1 };
	}
	else
	{
		roiLeft = { real(leftEndPoint), imag(leftEndPoint) };
		roiRight = { real(rightEndPoint), imag(rightEndPoint) };
	}
	
	//--------l------r-------*--------------------------		real(roi) = 1-real(rightEndPoint)..
	//--------l------r-------*-------R------L-----------						1-real(leftEndPoint)
	//
	//---------------l-------*--r-----------------------		real(roi) = real(rightEndPoint)..
	//---------------l----R--*--r----L------------------						1-real(leftEndPoint)
	//
	//---------------------l-*-------r------------------		real(roi) = 1-real(leftEndPoint)..
	//---------------R-----l-*-L-----r------------------						real(rightEndPoint)
	//
	//-----------------------*-------l---r--------------		real(poi) = real(leftEndPoint)..
	//-----------R---L-------*-------l---r--------------						real(rightEndPoint)

	// [roiLeft,roiRight] ready and valid now
	
	const arg_t realStep = { real(step), 0 };
	const arg_t imagStep = { 0, imag(step) };

	int64_t computedValuesSoFar = 0;
	
	for (auto z = leftEndPoint; imag(z) <= imag(rightEndPoint); z += imagStep)
	{
		for (z = { real(leftEndPoint), imag(z) }; real(z) <= real(rightEndPoint); z += realStep)
		{
			if (USE_PRECISION_LEVELSETS)
			{
				// under construction
			}
			else
			{
				auto CVZitrs = RiemannZetaFnCVZitrsNeeded(z, decimal_precision);
				auto DEFitrs = RiemannZetaFnDEFitrsNeeded(z, decimal_precision);
				auto zeta_z = CVZitrs < DEFitrs ? CalculateRiemannZetaFnCVZ(z,CVZitrs) : CalculateRiemannZetaFnDEF(z, DEFitrs);
				
				if (direct_methods_enough)
				{
					// if whole region in { Re z, Im z >= 0 } --> take no care of reflections
					valuepair_ready(z,zeta_z);
					++computedValuesSoFar;
				}
				else if (symmetric_to_center_both)
				{
					// if whole region is symmetric to the reflection point --> generate all reflections
					arg_t refl_z, refl_zeta_z;
					arg_t conj_z, conj_zeta_z;
					arg_t conj_refl_z, conj_refl_zeta_z;
					CalculateRiemannZetaFnReflectionFormula(z, zeta_z, refl_z, refl_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(z, zeta_z, conj_z, conj_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(refl_z, refl_zeta_z, conj_refl_z, conj_refl_zeta_z);
					valuepair_ready(z,zeta_z);
					valuepair_ready(refl_z, refl_zeta_z);
					valuepair_ready(conj_z, conj_zeta_z);
					valuepair_ready(conj_refl_z, conj_refl_zeta_z);
					computedValuesSoFar += 4;
				}
				else if (symmetric_to_center_real)
				{
					// horizontally symmetric case --> add both original and reflection in horizontal direction, but check for boundaries in vertical
					arg_t refl_z, refl_zeta_z;
					arg_t conj_z, conj_zeta_z;
					arg_t conj_refl_z, conj_refl_zeta_z;
					CalculateRiemannZetaFnReflectionFormula(z, zeta_z, refl_z, refl_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(z, zeta_z, conj_z, conj_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(refl_z, refl_zeta_z, conj_refl_z, conj_refl_zeta_z);
					if (imag(leftEndPoint) <= imag(z) && imag(z) <= imag(rightEndPoint))
					{
						valuepair_ready(z,zeta_z);
						valuepair_ready(refl_z, refl_zeta_z);
						computedValuesSoFar += 2;
					}
					if (imag(leftEndPoint) <= imag(conj_z) && imag(conj_z) <= imag(rightEndPoint))
					{
						valuepair_ready(conj_z, conj_zeta_z);
						valuepair_ready(conj_refl_z, conj_refl_zeta_z);
						computedValuesSoFar += 2;
					}
				}
				else if (symmetric_to_center_imag)
				{
					// vertically symmetric case --> add both original and conjugate if needed, but check for boundaries in horizontal direction
					arg_t refl_z, refl_zeta_z;
					arg_t conj_z, conj_zeta_z;
					arg_t conj_refl_z, conj_refl_zeta_z;
					CalculateRiemannZetaFnReflectionFormula(z, zeta_z, refl_z, refl_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(z, zeta_z, conj_z, conj_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(refl_z, refl_zeta_z, conj_refl_z, conj_refl_zeta_z);
					if (real(leftEndPoint) <= real(z) && real(z) <= real(rightEndPoint))
					{
						valuepair_ready(z,zeta_z);
						valuepair_ready(conj_z, conj_zeta_z);
						computedValuesSoFar += 2;
					}
					if (real(leftEndPoint) <= real(refl_z) && real(refl_z) <= real(rightEndPoint))
					{
						valuepair_ready(refl_z, refl_zeta_z);
						valuepair_ready(conj_refl_z, conj_refl_zeta_z);
						computedValuesSoFar += 2;
					}
				}
				else
				{
					// general case --> calculate all four values, but add only those requested originally
					arg_t refl_z, refl_zeta_z;
					arg_t conj_z, conj_zeta_z;
					arg_t conj_refl_z, conj_refl_zeta_z;
					auto ptInRegion = [&leftEndPoint, &rightEndPoint](arg_t point) {
						return real(leftEndPoint) <= real(point) && real(point) <= real(rightEndPoint) &&
							imag(leftEndPoint) <= imag(point) && imag(point) <= imag(rightEndPoint);
					};

					CalculateRiemannZetaFnReflectionFormula(z, zeta_z, refl_z, refl_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(z, zeta_z, conj_z, conj_zeta_z);
					CalculateRiemannZetaFnConjugateFormula(refl_z, refl_zeta_z, conj_refl_z, conj_refl_zeta_z);

					if (ptInRegion(z))
					{
						valuepair_ready(z, zeta_z);
						++computedValuesSoFar;
					}
					if (ptInRegion(refl_z))
					{
						valuepair_ready(refl_z, refl_zeta_z);
						++computedValuesSoFar;
					}
					if (ptInRegion(conj_z))
					{
						valuepair_ready(conj_z, conj_zeta_z);
						++computedValuesSoFar;
					}
					if (ptInRegion(conj_refl_z))
					{
						valuepair_ready(conj_refl_z, conj_refl_zeta_z);
						++computedValuesSoFar;
					}
				}					
			} // if (!USE_PRECISION_LEVELSETS)
			if (progress)
				progress(computedValuesSoFar);
		} // for real(z)
	} // for imag(z)
}


} // namespace numericalzeta