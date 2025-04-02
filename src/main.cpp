#include "../include/IPM.hpp"
#include "../include/cube_samp.hpp"
#include <cstdlib>   // for std::getenv
#include <string>

std::string get_home_directory() {
    const char* home = std::getenv("HOME");
    return home ? std::string(home) : "";
}


int main() 
{

    int n = 3;
    int N = 2000;
    double step = 0.025;

    Eigen::VectorXd a; // lower bounds
    Eigen::VectorXd b; // upper bounds
    Eigen::VectorXd c; // origin
    // Eigen::MatrixXd X; // points generated X(N, n)

    init_cube(n, a, b, c);
    auto [dist, X] = sample_n_cube_sparse(n, N, step, a, b, c);
    auto [dist_sparse, X_sparse] = sample_n_cube_sparse(n, N, step, a, b, c);

    std::cout << "expected dist. from origin - plain samp.: " << dist << "\n";
    std::cout << "expected dist. from origin - sparse samp.: " << dist_sparse << "\n";
    
    write_samples_to_file(X_sparse, "../points_sparse.txt");
    write_samples_to_file(X, "../points.txt");

    write_samples_to_file(X_sparse, get_home_directory()+"/Documents/Uni/comb-opt/set-cover/points_sparse.txt");
    write_samples_to_file(X, get_home_directory()+"/Documents/Uni/comb-opt/set-cover/points.txt");

    // HPolytopeType cube = init_cube_volesti(n, a, b, c);

    // std::cout << cube.get_vec() << "\n";
    // std::cout << cube.get_mat() << "\n";


    // LP LP_data;
    // State state;
    // Parameters params;

	// read_problem("../LP.txt", LP_data);
    // std::cout << "A\n" << LP_data.A << std::endl;
    // std::cout << "\nb\n" << LP_data.b << std::endl;
    // std::cout << "\nc\n" << LP_data.c << std::endl;

	// Primal_Dual_MPC(LP_data, state, params);
    // std::cout << "\nx\n" << state.x << std::endl;


	return 0;
}
