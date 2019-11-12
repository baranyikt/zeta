# Zeta
## Purpose
Goal of this is to calculate as many values of the [Riemann Zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function) as precise and as fast as possible -- in this order. It was made with two priciples in mind: one that it should be precise -- mathematically and by means of arbitrary numberical precision -- , two that it should be able to compute as many values as possible from given amount of work: recycling calculations, saving partial results and using reflection formulas to produce as much results as it can. In other words, it is optimized for creating a big image or table of Zeta values.

## Example usage

##### Most basic example, populate a large table of Zeta values in 0.3+0i..0.6+100i with stepping 0.0001+0.0001i (default) using default decimal precision (std::numeric_limits<arg_base_t>::digits10) and no progress callback:
```c++
#include <zeta.h>
using namespace numericalzeta;
std::vector<std::vector<arg_t>> resultMtx;
value_ready_func_t populatorFunction = ValueReadyFuncFactoryWithResultMatrix({0.3,0.0}, {0.6,100.0}, resultMtx);
CalculateRiemannZetaFn({0.3,0.0}, {0.6,100.0}, populatorFunction);
```
ValueReadyFuncFactoryWithResultMatrix() will call resultMtx.resize() to preinitialize resultMtx as two-dimensional vector to hold as many as 
$$[(real(rightEndPoint) - real(leftEndPoint)) / real(step)) + 1] * [(imag(rightEndPoint) - imag(leftEndPoint)) / imag(step)) + 1]$$
arg_t values. (Note that in this form it'll need ~24Gb of free memory for doing that.)

##### Same, but overriding defaults with chosen values:
```c++
std::vector<arg_t> resultMtx2;
const arg_t leftEndPoint	= {	0.3,	0.0		};
const arg_t rightEndPoint	= {	0.6,	100.0	};
const arg_t step			= { 0.01,	0.01 	};
uint32_t requestedDigits	= 7;
value_ready_func_t populatorFunction2 = ValueReadyFuncFactoryWithResultMatrix(leftEndPoint, rightEndPoint, resultMtx2, step);
auto progressFunction = [](int64_t computedSoFar) { DisplayPercent(static_cast<float>(computedSoFar)/resultMtx2.size() * 100); }
CalculateRiemannZetaFn(leftEndPoint, rightEndPoint, populatorFunction2, step, requestedDigits, progressFunction);
```
Note: progress function is given the number of values calculated so far, we need to divide that by the total number of values to get the percentage of the work done.

##### You can use your own function to collect calculated values, of course:
```c++
auto myValueCollector = [](const arg_t& x, const arg_t& y) { 
	arg_base_t hue, saturation, lightness;
	ConvertComplexToHSL(y, hue, saturation, lightness);
	SetPixelColor(scale(real(x)),scale(imag(x)),HSLtoRGB(hue, saturation, lightness));
}
CalculateRiemannZetaFn({0.3,0.0}, {0.6,100.0}, myValueCollector);
```
##### Some functions doing internal calculations are available, too:
```c++
// Get zeta(0.6+100i) using Cohen-Rodriguez Villegas-Zagier method with at least 8 digits precision
arg_t zetaResultCVZ = CalculateRiemannZetaFnCVZ({0.6,100.0}, RiemannZetaFnCVZitrsNeeded({0.6,100.0},8));

// Get zeta(0.3+100i) using Zeta series definition with at least 7 digits precision
arg_t zetaResultDEF = CalculateRiemannZetaFnDEF({0.3,100.0}, RiemannZetaFnCVZitrsNeeded({0.3,100.0},7));
```

## TODO
- implement speed of convergence function of Zeta series definition method for Im z != 0
- arbitrary precision gamma function for reflection formulas
- implement USE_PRECISION_LEVELSETS -- hardcoded precision boundaries for different methods in order to be able to choose the best method without the need of computing the number of required iterations for given precision at every point
- add Odlyzko--Schönhage method
- user-defined literals

## Links & references
[Xavier Gourdon, Pascal Sebah: Numerical evaluation of the Riemann Zeta-function (2003)](http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.pdf)
[A. Odlyzko, A. Schönhage: Fast algorithms for multiple evaluations of the Riemann Zeta-function (1988)](http://www.dtc.umn.edu/~odlyzko/doc/arch/fast.zeta.eval.pdf)
[H. Cohen, F. Rodriguez Villegas, D. Zagier: Convergence acceleration of alternating series (1991)](https://people.mpim-bonn.mpg.de/zagier/files/exp-math-9/fulltext.pdf)
[J. M. Borwein, D. M. Bradley, R. E. Crandall: Computational strategies for the Riemann zeta function](https://cr.yp.to/bib/2000/borwein.pdf)
[Riemann Zeta function](https://en.wikipedia.org/wiki/Riemann_zeta_function)
[Lanczos approximation for fixed precision gamma function](https://en.wikipedia.org/wiki/Lanczos_approximation)
[G. R. Pugh: An Analysis of the Lanczos Gamma Approximation (PhD thesis)](https://web.viu.ca/pughg/phdThesis/phdThesis.pdf)

##### Monographs on the topic:
Aleksandar Ivic: The Riemann Zeta-Function (1985, 2003)
E. C. Titchmarsh: The theory of the Riemann Zeta-Function (1951, 1986)
H. M. Edwards: Riemann's Zeta Function (1974)

