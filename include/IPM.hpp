#ifndef IPM_HPP
#define IPM_HPP

#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

#define TOL 1e-7
#define MAX_ITERS 3

struct State{
    Eigen::VectorXd x, y, s;
    Eigen::MatrixXd X, S, M, X_inv, S_inv, M_inv;
    Eigen::VectorXd rd, rp, rc, rhs;
    Eigen::VectorXd Dx_p, Dy_p, Ds_p, Dx, Dy, Ds;
    Eigen::VectorXd e;
};

struct Parameters{
    double mu, alpha_p_p, alpha_d_p, alpha_p, alpha_d, eta, sigma;
};

struct LP{
    Eigen::MatrixXd A, A_T;
    Eigen::VectorXd b, c;
    int n, m;
};

void read_problem (std::string filepath, LP &LP_data);

void initialize_algorithm (LP LP_data, State &state);

void predictor_step (LP LP_data, State &state, Parameters &params);

void corrector_step (LP LP_data, State &state, Parameters &params);

void Primal_Dual_MPC (LP LP_data, State &state, Parameters &params);

#endif //IPM_HPP