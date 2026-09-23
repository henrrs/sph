#include "sph_kernels.h"

#include "parameters.h"

namespace {

using namespace sim;

[[nodiscard]] inline double square(double value) noexcept {
	return value * value;
}

[[nodiscard]] inline double cube(double value) noexcept {
	return value * value * value;
}

[[nodiscard]] inline double fourthPower(double value) noexcept {
	const double value2 = square(value);
	return value2 * value2;
}

[[nodiscard]] inline double fifthPower(double value) noexcept {
	return fourthPower(value) * value;
}

[[nodiscard]] inline double smoothingLengthCubed() noexcept {
	return h * h * h;
}

} // namespace

double euclideanDistance(const Vec3& lhs, const Vec3& rhs) {
	return norm(rhs - lhs);
}

double cubicSplineKernel(double r) {
	const double q = r / sim::h;
	double wij = 0.0;
	const double alphaValue = 1.0 / (sim::pi * smoothingLengthCubed());

	if ((q >= 0.0) && (q <= 1.0)) {
		wij = 1.0 - (1.5 * square(q) * (1.0 - (q / 2.0)));
	}
	else if ((q > 1.0) && (q <= 2.0)) {
		wij = 0.25 * cube(2.0 - q);
	}

	return alphaValue * wij;
}

double kernelFunction(double r) {
	const double q = r / sim::h;
	double wij = 0.0;
	const double alphaValue = 3.0 / (359.0 * sim::pi * smoothingLengthCubed());

	if ((q >= 0.0) && (q <= 1.0)) {
		wij = fifthPower(3.0 - q) - (6.0 * fifthPower(2.0 - q)) + (15.0 * fifthPower(1.0 - q));
	}
	else if ((q > 1.0) && (q <= 2.0)) {
		wij = fifthPower(3.0 - q) - (6.0 * fifthPower(2.0 - q));
	}
	else if ((q > 2.0) && (q <= 3.0)) {
		wij = fifthPower(3.0 - q);
	}

	return alphaValue * wij;
}

Vec3 gradKernelFunction(double r, const Vec3& xij) {
	const double alphaValue = 3.0 / (359.0 * sim::pi * smoothingLengthCubed());
	Vec3 grad = { 0.0, 0.0, 0.0 };
	const double q = r / sim::h;
	double dfdq = 0.0;

	if (r <= 0.00001 || q >= 3.0) {
		return grad;
	}

	if (q < 1.0) {
		dfdq = (-5.0 * fourthPower(3.0 - q)) + (30.0 * fourthPower(2.0 - q)) - (75.0 * fourthPower(1.0 - q));
	}
	else if (q < 2.0) {
		dfdq = (-5.0 * fourthPower(3.0 - q)) + (30.0 * fourthPower(2.0 - q));
	}
	else {
		dfdq = -5.0 * fourthPower(3.0 - q);
	}

	const double gradientScale = alphaValue * dfdq / (sim::h * r);
	return gradientScale * xij;
}
