// vector2d.hpp
#ifndef VECTOR2D_HPP
#define VECTOR2D_HPP

#include <cmath>

struct Vector2D {
    double x, y;

    Vector2D() : x(0.0), y(0.0) {}
    Vector2D(double x, double y) : x(x), y(y) {}

    Vector2D operator-(const Vector2D& other) const { return {x - other.x, y - other.y}; }
    Vector2D& operator-=(const Vector2D& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    Vector2D operator+(const Vector2D& other) const { return {x + other.x, y + other.y}; }
    Vector2D& operator+=(const Vector2D& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vector2D operator/(const double k) const { return {x / k, y / k}; }
    Vector2D operator*(const double k) const { return {x * k, y * k}; }

    // Friend function for k * v
    friend Vector2D operator*(const double k, const Vector2D& v) { return v * k; }
};

inline double magnitude(const Vector2D& v) { return std::sqrt(v.x * v.x + v.y * v.y); }

inline double distance(const Vector2D& a, const Vector2D& b) { return magnitude(a - b); }

#endif  // VECTOR2D_HPP