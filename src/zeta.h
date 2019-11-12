#include <complex>
#include <functional>
#include <vector>
#include <limits>

namespace numericalzeta {

using arg_base_t = double;
using intrmed_base_t = long double;
using arg_t = std::complex<arg_base_t>;
using intrmed_t = std::complex<intrmed_base_t>;
using value_ready_func_t = std::function<void(const arg_t&,const arg_t&)>;
using progress_func_t = std::function<void(int64_t)>;
using intrmed_int_t = int64_t;

const arg_t DefaultStep = { 0.0001, 0.0001 };

constexpr bool USE_PRECISION_LEVELSETS = false;
constexpr double FLOATINGPOINT_EQ_TOLERANCE = 10e-7;

value_ready_func_t ValueReadyFuncFactoryWithResultMatrix(arg_t leftEndPoint, arg_t rightEndPoint, std::vector<std::vector<arg_t>>& mtx, arg_t step = DefaultStep);

intrmed_int_t RiemannZetaFnCVZitrsNeeded(arg_t z, intrmed_int_t decimal_precision);
arg_t CalculateRiemannZetaFnCVZ(arg_t z, intrmed_int_t iterations);

intrmed_int_t RiemannZetaFnDEFitrsNeeded(arg_t z, intrmed_int_t decimal_precision);
arg_t CalculateRiemannZetaFnDEF(arg_t z, intrmed_int_t iterations);

void CalculateRiemannZetaFn(arg_t leftEndPoint, arg_t rightEndPoint, value_ready_func_t valuepair_ready,
	arg_t step = DefaultStep, uint32_t decimal_precision = std::numeric_limits<arg_base_t>::digits10, progress_func_t progress = {});

} // namespace numericalzeta