#include "cube_samp.hpp"

void init_cube (int n, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c) 
{
    std:: cout << "Cube dimension: " << n << std::endl;
    a.resize(n);
    b.resize(n);

    for (int i=0; i<n; i++)
    {
        a(i) = -1.0;
        b(i) = 1.0;
    }

    // calc. origin
    c = (a + b) / 2.0;
}


/* using this for sample_n_cube function */
double expected_dist_from_origin (int N, Eigen::MatrixXd X, Eigen::VectorXd c)
{
    double dist = 0.0;

    for (int i = 0; i<N; i++)
    {
        Eigen::VectorXd x = X.row(i);
        dist += (x - c).norm(); // euclidean 
    }


    return dist / N;
}


std::pair<double, Eigen::MatrixXd> sample_n_cube (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd c)
{   
    Eigen::MatrixXd X(N, n);
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> index_dist(0, n - 1);
    std::uniform_int_distribution<> bound_dist(0, 1);

    std::set<std::string> seen;

    for (int i=0; i<N; i++)
    {
        Eigen::VectorXd x(n);

        // fix coordinate
        int fixed_coordinate = index_dist(gen);

        // fix side
        bool upper = bound_dist(gen);

        x(fixed_coordinate) = upper ? b(fixed_coordinate) : a(fixed_coordinate);

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
        }

        std::string key = vector_to_string(x);
        if (seen.find(key) == seen.end())
        {
            X.row(i) = x;
            // std::cout << x.transpose() << "\n";
            // std::cout << i << "\n";
            seen.insert(key);
        }
        else
            i--;
    }

    double expected_dist = expected_dist_from_origin(N, X, c);
    std::pair<double, Eigen::MatrixXd> p = std::make_pair(expected_dist, X);

    return p;
}



std::pair<double, std::vector<Eigen::SparseVector<double>>> sample_n_cube_sparse (int n, int N, double step, Eigen::VectorXd a, Eigen::VectorXd b, Eigen::VectorXd c)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_int_distribution<> coord_dist(0, n - 1);
    std::uniform_int_distribution<> bound_dist(0, 1);

    std::set<std::string> seen;
    std::vector<Eigen::SparseVector<double>> samples;
    samples.reserve(N);

    double total_distance = 0.0;

    for (int i=0; i<N; i++)
    {
        Eigen::SparseVector<double> x(n);

        int fixed_coordinate = coord_dist(gen);
        bool upper = bound_dist(gen);
        x.insert(fixed_coordinate) = upper ? b(fixed_coordinate) : a(fixed_coordinate);

        for (int j = 0; j < n; ++j) {
            if (j == fixed_coordinate) continue;

            std::vector<double> grid;
            double from = a(j);
            while (from <= b(j) - 1e-10) {
                grid.push_back(from);
                from += step;
            }

            std::uniform_int_distribution<> grid_index(0, int(grid.size()) - 1);
            x.insert(j) = grid[grid_index(gen)];
        }

        std::string key = vector_to_string(x);
        if (seen.find(key) == seen.end()) {
            seen.insert(key);
            samples.push_back(x);

            double dist_sq = 0.0;
            for (Eigen::SparseVector<double>::InnerIterator it(x); it; ++it) {
                double d = it.value() - c(it.index());
                dist_sq += d * d;
            }
            total_distance += std::sqrt(dist_sq);
        }
        else
            i--;
    }

    double expected_dist = total_distance / N;
    std::pair<double, std::vector<Eigen::SparseVector<double>>> p = std::make_pair(expected_dist, samples);

    return p;
}



HPolytopeType init_cube_volesti (int n, Eigen::VectorXd& a, Eigen::VectorXd& b, Eigen::VectorXd& c)
{   
    // define lower and upper bounds
    a.resize(n);
    b.resize(n);
    for (int i=0; i<n; i++)
    {
        a(i) = -1;
        b(i) = 1;
    }

    c = (a + b) / 2.0;

    Eigen::MatrixXd A(2*n, n);
    A.topRows(n) = Eigen::MatrixXd::Identity(n, n);
    A.bottomRows(n) = -Eigen::MatrixXd::Identity(n, n);

    Eigen::VectorXd b_(2*n);
    b_.head(n) = b;
    b_.tail(n) = -a;

    return HPolytopeType(n, A, b_);
}



void volesti_samp (int n, int N, int walk_len, HPolytopeType cube, Eigen::VectorXd starting_point)
{
    typedef BoostRandomNumberGenerator<boost::mt19937, double> RNGType;
    RNGType rng(n);

    std::vector<Point> samples;

    /* add this as a parameter to the function */
    Point p(n), p1, p2;
    for (int i=0; i<n; ++i) p.set_coord(i, starting_point(i));
    // p.print();

    BRDHRWalk::Walk<HPolytopeType, RNGType> walk(cube, p, rng);

    for (int i=0; i<N; i++)
    {
        walk.apply(cube, p1, p2, walk_len, rng);
        samples.push_back(p2);
        p2.print();
    }
}