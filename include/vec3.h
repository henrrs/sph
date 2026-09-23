#pragma once

#include <cstddef>
#include <cmath>

struct Vec3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;

	constexpr double operator[](std::size_t index) const noexcept {
		return index == 0 ? x : (index == 1 ? y : z);
	}

	constexpr double& operator[](std::size_t index) noexcept {
		return index == 0 ? x : (index == 1 ? y : z);
	}
};

[[nodiscard]] inline constexpr Vec3 operator+(const Vec3& lhs, const Vec3& rhs) noexcept {
	return { lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
}

[[nodiscard]] inline constexpr Vec3 operator-(const Vec3& lhs, const Vec3& rhs) noexcept {
	return { lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
}

[[nodiscard]] inline constexpr Vec3 operator*(const Vec3& value, double scale) noexcept {
	return { value.x * scale, value.y * scale, value.z * scale };
}

[[nodiscard]] inline constexpr Vec3 operator*(double scale, const Vec3& value) noexcept {
	return value * scale;
}

inline constexpr Vec3& operator+=(Vec3& lhs, const Vec3& rhs) noexcept {
	lhs.x += rhs.x;
	lhs.y += rhs.y;
	lhs.z += rhs.z;
	return lhs;
}

inline constexpr Vec3& operator-=(Vec3& lhs, const Vec3& rhs) noexcept {
	lhs.x -= rhs.x;
	lhs.y -= rhs.y;
	lhs.z -= rhs.z;
	return lhs;
}

[[nodiscard]] inline constexpr double dot(const Vec3& lhs, const Vec3& rhs) noexcept {
	return (lhs.x * rhs.x) + (lhs.y * rhs.y) + (lhs.z * rhs.z);
}

[[nodiscard]] inline double norm(const Vec3& value) noexcept {
	return std::sqrt(dot(value, value));
}
