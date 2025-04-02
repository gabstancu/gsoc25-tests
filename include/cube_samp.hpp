#ifndef CUBE_SAMP_HPP
#define CUBE_SAMP_HPP

#include "Eigen/Eigen"
#include "cartesian_geom/cartesian_kernel.h"
#include "convex_bodies/hpolytope.h"
#include "generators/known_polytope_generators.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>
#include <set>

typedef Cartesian<double> Kernel;
typedef Kernel::Point Point;
typedef HPolytope<Point> HPolytopeType;


template <typename VectorType>
std::string vector_to_string(const VectorType& x) {
    std::ostringstream oss;
    oss.precision(10);

    if constexpr (std::is_same_v<VectorType, Eigen::SparseVector<double>>) {
        for (typename Eigen::SparseVector<double>::InnerIterator it(x); it; ++it) {
            oss << it.index() << ":" << it.value() << " ";
        }
    } else {
        for (int i = 0; i < x.size(); ++i) {
            oss << x(i);
            if (i != x.size() - 1) oss << " ";
        }
    }

    return oss.str();
}

template <typename>
struct always_false : std::false_type {};

template <typename T>
void print_samples(const T& samples) {
    
    if constexpr (std::is_same_v<T, Eigen::MatrixXd>) {
        for (int i = 0; i < samples.rows(); ++i) {
            std::cout << "Sample " << i << ": ";
            for (int j = 0; j < samples.cols(); ++j) {
                std::cout << samples(i, j);
                if (j < samples.cols() - 1) std::cout << " ";
            }
            std::cout << "\n";
        }
    }
    else if constexpr (std::is_same_v<T, std::vector<Eigen::SparseVector<double>>>) {
        for (size_t i = 0; i < samples.size(); ++i) {
            std::cout << "Sample " << i << ": ";
            for (Eigen::SparseVector<double>::InnerIterator it(samples[i]); it; ++it) {
                std::cout << it.value() << " ";
            }
            std::cout << "\n";
        }
    }
    else {
        static_assert(always_false<T>::value, "Unsupported type passed to print_samples!");
    }
}


template <typename SampleType>
void write_samples_to_file(const SampleType& samples, const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    if constexpr (std::is_same_v<SampleType, Eigen::MatrixXd>) {
        for (int i = 0; i < samples.rows(); ++i) {
            for (int j = 0; j < samples.cols(); ++j) {
                file << samples(i, j);
                if (j < samples.cols() - 1) file << " ";
            }
            file << "\n";
        }
    }
    else if constexpr (std::is_same_v<SampleType, std::vector<Eigen::SparseVector<double>>>) {
        for (const auto& x : samples) {
            Eigen::VectorXd dense = Eigen::VectorXd(x);  // convert to full vector
            for (int j = 0; j < dense.size(); ++j) {
                file << dense(j);
                if (j < dense.size() - 1) file << " ";
            }
            file << "\n";
        }
    }
    else {
        static_assert(always_false<SampleType>::value, "Unsupported type passed to write_samples_to_file!");
    }

    file.close();
    std::cout << "Wrote samples to " << filename << "\n";
}


void init_cube (int n, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c);


double expected_dist_from_origin (int N, Eigen::MatrixXd X, Eigen::VectorXd c);


std::pair<double, Eigen::MatrixXd> sample_n_cube (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd c);


std::pair<double, std::vector<Eigen::SparseVector<double>>> sample_n_cube_sparse (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd c);


HPolytopeType init_cube_volesti (int n, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c);


void volesti_samp (int n, int N, double step, const Eigen::VectorXd& a, const Eigen::VectorXd& b);


#endif // CUBE_SAMP_HPP