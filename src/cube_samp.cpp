#include "cube_samp.hpp"


// initializes a random n-dimensional cube: (-1, 1)^n by default
// initializes lower and upper bounds, a and b
void init_cube (int n, Eigen::VectorXd& a, Eigen::VectorXd& b) 
{
    std:: cout << "Cube dimension: " << n << std::endl;
    a.resize(n);
    b.resize(n);

    for (int i=0; i<n; i++)
    {
        a(i) = -1.0;
        b(i) = 1.0;
    }
}

double expected_dist_from_origin (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b,  Eigen::MatrixXd X)
{
    double dist = 0.0;



    return dist;
}

// a: lower bounds
// b: upper bounds
void sample_n_cube (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b,  Eigen::MatrixXd& X)
{   

    X.resize(N, n);

    for (int i=0; i<N; i++)
    {
        Eigen::VectorXd x(n);

        std::random_device rd;
        std::mt19937 gen(rd());

        // fix coordinate
        std::uniform_int_distribution<> index_dist(0, n - 1);
        int fixed_coordinate = index_dist(gen);

        // fix side
        std::uniform_int_distribution<> bound_dist(0, 2);
        bool upper = bound_dist(gen);

        for (int j = 0; j<n; j++)
        {   
            if (j==fixed_coordinate) continue;

            /* generate grid according to step */
            std::vector<double> grid;
            double from = a(j);

            while (from <= b(j) - 1e-10)
            {
                grid.push_back(from);
                // std:: cout << from << " ";
                from += step;
            }

            // pick randomly from grid
            std::uniform_int_distribution<> grid_index(0, int(grid.size()) - 1);
            int grid_idx = grid_index(gen);
            x[j] = grid[grid_idx];
            X.row(i) = x;
            
        }
    }
}