#include "IPM.hpp"

void read_problem(const std::string filepath, LP &LP_data)
{
    
    std::fstream data(filepath);
    std::vector<std::string> lines;
    std::string line;

    if(!data)
    {
        std::cerr << "Error opening " << filepath << std::endl;
        exit(1);
    }

    while(std::getline(data, line))
    {   
        lines.push_back(line);
    }

    // get num_vars and num_ctrs
    std::string num;
    std::vector<int> v;
    for(int i=0; i<lines[0].size(); i++)
    {   
        if(lines[0][i]==' ')
        {   
            v.push_back(std::stoi(num));
            num = "";
            continue;
        }
        num+=lines[0][i];
    }

    num = "";

    int num_vars = v[0];
    int num_ctrs = v[1];
    LP_data.A.resize(num_ctrs, num_vars);
    LP_data.b.resize(num_ctrs, 1);
    LP_data.c.resize(num_vars, 1);

    int index = 0;
    // coefficients line
    for(int i=0; i<lines[1].size(); i++)
    {
        if(lines[1][i]==' ')
        {   
            LP_data.c(index++, 0) = std::stoi(num);
            num = "";
            continue;
        }
        num+=lines[1][i];
    }
    
    // constraints lines
    std::vector<int> m;
    for(int i = 2; i<lines.size(); i++)
    {
        for(int j=0; j<lines[i].size(); j++)
        {  
            if(lines[i][j]==' ')
            {   
                m.push_back(std::stoi(num));
                num = "";
                continue;
            }
            num+=lines[i][j];

        }
    }


    int step = 0;
    int b_index = 0;
    int a_row = 0;
    int a_col = 0;
    for(int i=0; i<m.size(); i++)
    {   
        if(step==num_ctrs)
        {
            step = 0;
            LP_data.b(b_index++, 0) = m[i];
            a_row++;
            a_col = 0;
            continue;
        }
        step++;
        LP_data.A(a_row, a_col++) = m[i];

    }

    LP_data.A_T = LP_data.A.transpose();

    data.close();
}


void initialize_algorithm (LP LP_data, State &state)
{
    /* comp. starting point using Mehrotra's heuristic */
    Eigen::MatrixXd H = LP_data.A * LP_data.A_T;
    Eigen::VectorXd x_ = LP_data.A_T * H.inverse() * LP_data.b;
    Eigen::VectorXd y_ = H.inverse() * LP_data.A * LP_data.c;
    Eigen::VectorXd s_ = LP_data.c - LP_data.A_T * y_;

    state.e = Eigen::VectorXd::Ones(x_.size());

    double dx = std::max(-1.5 * x_.minCoeff(), 0.0);
    double ds = std::max(-1.5 * s_.minCoeff(), 0.0);
    double dx_ = dx + 0.5 * ( ( (x_ + dx * state.e).transpose() * ( s_ + ds * state.e ) ) / (s_.sum() + dx) ).coeff(0);
    double ds_ = ds + 0.5 * ( ( (x_ + dx * state.e).transpose() * ( s_ + ds * state.e ) ) / (x_.sum() + ds) ).coeff(0);

    state.x = x_ + dx_ * state.e;
    state.y = y_;
    state.s = s_ + ds_ * state.e;

    std::cout << "\nstarting point:" << std::endl;
    std::cout << "\nx\n" << state.x << std::endl; 
    std::cout << "\ny\n" << state.y << std::endl; 
    std::cout << "\ns\n" << state.s << std::endl; 
}


void predictor_step (LP LP_data, State &state, Parameters &params)
{
    std::cout << "*************** predictor step ***************" << std::endl;
    state.rhs = state.rp + LP_data.A * state.S_inv * (-state.rc + state.X * state.rd);
    std::cout << "\nrhs:\n" << state.rhs << std::endl; 

    state.Dy_p = state.M_inv * state.rhs;
    state.Ds_p = state.rd - LP_data.A_T * state.Dy_p;
    state.Dx_p = state.S_inv * (state.rc - state.X * state.Ds_p);

    std::cout << "\nDy_p:\n" << state.Dy_p << std::endl; 
    std::cout << "Ds_p:\n" << state.Ds_p << std::endl;
    std::cout << "Dx_p:\n" << state.Dx_p << std::endl;


    // calculate alpha for both the primal and the dual problem
    std::vector<double> x_Dx = {};
    for (int j=0; j<state.Dx_p.size(); j++)
    {   
        if(state.Dx_p(j)>0) x_Dx.push_back(state.x(j)/state.Dx_p(j));
    }
    
    if (!x_Dx.empty())
        params.alpha_p_p = std::min(1.0, *std::min_element(x_Dx.begin(), x_Dx.end()));
    else
        params.alpha_p_p = 1.0;
    
    std::vector<double> s_Ds = {};
    for (int j=0; j<state.Ds_p.size(); j++)
    {   
        if(state.Ds_p(j)>0) s_Ds.push_back(state.s(j)/state.Ds_p(j));
    }
    
    if (!s_Ds.empty())
        params.alpha_d_p = std::min(1.0, *std::min_element(s_Ds.begin(), s_Ds.end()));
    else 
        params.alpha_d_p = 1.0;
    
    std::cout << "\nalpha_primal_p: " << params.alpha_p_p << " alpha_dual_p: " << params.alpha_d_p << std::endl;
}


void corrector_step (LP LP_data, State &state, Parameters &params)
{   
    std::cout << "*************** corrector step ***************" << std::endl;
    state.e = Eigen::VectorXd::Ones(state.rc.size());

    state.rc = state.rc - params.sigma * params.mu * state.e + state.Dx_p.asDiagonal() * state.Ds_p;
    state.rhs = state.rp + LP_data.A * state.S_inv * (-state.rc + state.X * state.rd);
    std::cout << "\nrc:\n" << state.rc << std::endl; 
    std::cout << "rhs:\n" << state.rhs << std::endl;

    state.Dy = state.M_inv * state.rhs;
    state.Ds = state.rd - LP_data.A_T * state.Dy;
    state.Dx = state.S_inv * (state.rc - state.X * state.Ds);

    std::cout << "\nDy:\n" << state.Dy << std::endl; 
    std::cout << "Ds:\n" << state.Ds << std::endl;
    std::cout << "Dx:\n" << state.Dx << std::endl;

    // calc. primal and dual steps
    params.eta = std::max(0.995, 1 - params.mu);

    std::vector<double> x_Dx = {};
    for (int j=0; j<state.Dx.size(); j++)
    {
        if (state.Dx(j)>0) x_Dx.push_back(state.x(j)/state.Dx(j));
    }

    if(!x_Dx.empty())
        params.alpha_p = std::min(1.0, params.eta * *std::min_element(x_Dx.begin(), x_Dx.end()));
    else
        params.alpha_p = 1.0;

        
    std::vector<double> s_Ds = {};
    for (int j=0; j<state.Ds.size(); j++)
    {
        if (state.Ds(j)>0) s_Ds.push_back(state.s(j)/state.Ds(j));
    }
    
    if(!s_Ds.empty())
        params.alpha_d = std::min(1.0, params.eta * *std::min_element(s_Ds.begin(), s_Ds.end()));
    else
        params.alpha_d = 1.0;
    
    std::cout << "\nalpha_primal: " << params.alpha_p << " alpha_dual: " << params.alpha_d << std::endl;
}


void Primal_Dual_MPC (LP LP_data, State &state, Parameters &params)
{   
    initialize_algorithm(LP_data, state);

    for (int i=0; i<MAX_ITERS; i++)
    {   
        std::cout << "------------------------------ iter. " << i+1 << " ------------------------------" << std::endl;
        /* set diagonals */
        state.X = state.x.asDiagonal();
        state.S = state.s.asDiagonal();
        state.S_inv = state.S.inverse();

        /* calc. residuals */
        state.rd = LP_data.A_T * state.y + state.s - LP_data.c;
        state.rp = LP_data.A * state.x - LP_data.b;
        state.rc = state.X * state.s;

        std::cout << "\nrd:\n" << state.rd << std::endl; 
        std::cout << "rp:\n" << state.rp << std::endl;
        std::cout << "rc:\n" << state.rc << std::endl;

        params.mu = ((state.x.transpose() * state.s) / LP_data.n).coeff(0);
        std::cout << "mu: " << params.mu << std::endl;

        /* optimality check */
        if(std::max({params.mu, state.rp.norm(), state.rd.norm()}) <= TOL)
        {
            std::cout << "LP is optimal." << std::endl;
            break;
        }

        /* calc. M for both systems */
        state.M = LP_data.A * state.X * state.S_inv * LP_data.A_T;
        state.M_inv = state.M.inverse();

        predictor_step(LP_data, state, params);

        /* calc. sigma */
        params.sigma = std::pow(((state.x - params.alpha_p_p * state.Dx_p).transpose() * (state.s - params.alpha_d_p * state.Ds_p) / (LP_data.n * params.mu) ).coeff(0), 3);
        // std::cout << "sigma: " << params.sigma << std::endl;

        corrector_step(LP_data, state, params);

        /* update (x, y, s) */
        state.x = state.x - params.alpha_p * state.Dx;
        state.y = state.y - params.alpha_d * state.Dy;
        state.s = state.s - params.alpha_d * state.Ds;
    }
}