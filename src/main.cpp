#include "../include/IPM.hpp"
#include "../include/cube_samp.hpp"

int main() 
{

    int n = 5;
    int N = 1000;
    double step = 0.25;

    Eigen::VectorXd a; // lower bounds
    Eigen::VectorXd b; // upper bounds
    Eigen::MatrixXd X; // points generated X(N, n)

    init_cube(n, a, b);
    sample_n_cube(n, N, step, a, b, X);

    std::cout << X << "\n";

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
