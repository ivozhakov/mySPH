// vector.hpp
#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <cmath>
#include <array>  // Для хранения компонент

template<int Dim>
struct Vector {
    std::array<double, Dim> data;  // x, y (для Dim=2); x, y, z (для Dim=3)

    Vector() { data.fill(0.0); }
    explicit Vector(double val) { data.fill(val); }
    Vector(std::initializer_list<double> init) {
        std::copy(init.begin(), init.end(), data.begin());
    }

    double& operator[](size_t i) { 
        return data[i]; 
    }

    const double& operator[](size_t i) const { 
        return data[i]; 
    }

    Vector operator-(const Vector& other) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) res[i] = data[i] - other[i];
        return res;
    }
    Vector& operator-=(const Vector& other) {
        for (int i = 0; i < Dim; ++i) data[i] -= other[i];
        return *this;
    }
    Vector operator+(const Vector& other) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) res[i] = data[i] + other[i];
        return res;
    }
    Vector& operator+=(const Vector& other) {
        for (int i = 0; i < Dim; ++i) data[i] += other[i];
        return *this;
    }
    Vector operator/(double k) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) res[i] = data[i] / k;
        return res;
    }
    Vector operator*(double k) const {
        Vector res;
        for (int i = 0; i < Dim; ++i) res[i] = data[i] * k;
        return res;
    }

    friend Vector operator*(double k, const Vector& v) { return v * k; }
};

template<int Dim>
double magnitude(const Vector<Dim>& v) {
    double sum = 0.0;
    for (double val : v.data) sum += val * val;
    return std::sqrt(sum);
}

template<int Dim>
double squared_norm(const Vector<Dim>& v) {
    double sum = 0.0;
    for (double val : v.data) sum += val * val;
    return sum;
}

template<int Dim>
double distance(const Vector<Dim>& a, const Vector<Dim>& b) {
    return magnitude(a - b);
}

template<int Dim>
double dot(const Vector<Dim>& a, const Vector<Dim>& b) {
    double sum = 0.0;
    for (int i = 0; i < Dim; ++i) sum += a[i] * b[i];
    return sum;
}

#endif  // VECTOR_HPP