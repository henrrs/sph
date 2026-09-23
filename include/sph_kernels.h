#pragma once

#include "vec3.h"

double euclideanDistance(const Vec3& lhs, const Vec3& rhs);
double cubicSplineKernel(double r);
double kernelFunction(double r);
Vec3 gradKernelFunction(double r, const Vec3& xij);
