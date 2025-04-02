#ifndef CUBE_SAMP_HPP
#define CUBE_SAMP_HPP

#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <random>


void init_cube (int n, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c);

double expected_dist_from_origin (int N, Eigen::MatrixXd X, Eigen::VectorXd c);

void sample_n_cube (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::MatrixXd& X);


#endif // CUBE_SAMP_HPP