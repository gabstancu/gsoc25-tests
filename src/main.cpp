#include "../include/IPM.hpp"

int main() 
{

    LP LP_data;
    State state;
    Parameters params;

	read_problem("/Users/gabrielastancu/Documents/gsoc25/test-project/LP.txt", LP_data);
    std::cout << "A\n" << LP_data.A << std::endl;
    std::cout << "b\n" << LP_data.b << std::endl;
    std::cout << "c\n" << LP_data.c << std::endl;

	Primal_Dual_MPC(LP_data, state, params);

    std::cout << "x\n" << state.x << std::endl;


	return 0;
}
