#include "../inc/Bezier.h"
#include <iostream>

// Function to compute the Kronecker product
matrix_t kroneckerProduct(const matrix_t& A, const matrix_t& B) {
    // ChatGPT generated
    matrix_t result(A.rows() * B.rows(), A.cols() * B.cols());

    for (int i = 0; i < A.rows(); ++i) {
        for (int j = 0; j < A.cols(); ++j) {
            result.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
        }
    }

    return result;
}

matrix_t Bezier::S_matrix(int order) {
  matrix_t S;
  S.resize(order,order+1);
  S.setZero();
  for (int i = 0; i < order; i++) {
    S(i,i) = -order;
    S(i,i+1) = order;
  }
  return S;
}

matrix_t Bezier::R_matrix(int order) {
  matrix_t R;
  R.resize(order+1,order);
  R.setZero();

  for (int i = 0; i < order; i++) {
    R(i,i) = ((float) order-i)/order;
    R(i+1,i) = ((float) i+1)/order;
  }
  return R;
}

matrix_t Bezier::H_matrix(int order) {
  matrix_t H_matrix(order+1, order+1);
  H_matrix = (1./tau)*S.transpose()*R.transpose();
  return H_matrix;
}


matrix_t Bezier::H_vec(const matrix_t &H, const int m, const int order, const int gamma, const int power_of_H) {
    // ChatGPT generated
    matrix_t K_mat = K_matrix(order + 1, power_of_H + 1);

    // Initialize H_tmp as an empty matrix
    matrix_t H_tmp(m * H.cols() * (power_of_H + 1), m * H.rows());
    H_tmp.setZero();

    int row_offset = 0;
    for (int i = 0; i <= power_of_H; ++i) {
        matrix_t H_transformed = H.transpose().pow(i);
        matrix_t kron_prod = kroneckerProduct(H_transformed, matrix_t::Identity(m, m));
        H_tmp.block(row_offset, 0, kron_prod.rows(), kron_prod.cols()) = kron_prod;
        row_offset += kron_prod.rows();
    }

    // Compute the final H_vec_
    matrix_t H_vec_ = kroneckerProduct(K_mat, matrix_t::Identity(m, m)) * H_tmp;

    return H_vec_;
}

matrix_t Bezier::Pi(const vector_t& c, int m, int gamma) {
    // Define the matrix C
    matrix_t C(2, 4);
    C << c(0),  c(0), -c(0), -c(0),
         c(1), -c(1),  c(1), -c(1);

    int n = m * gamma;

    // Create the identity matrices A and B
    matrix_t A = matrix_t::Identity(n, n);
    matrix_t B = matrix_t::Identity(m, m);

    // Create Pi_1 and Pi_2 using Kronecker product
    matrix_t Pi_1 = kroneckerProduct(A, matrix_t::Ones(m, 1)).eval();
    matrix_t Pi_2 = kroneckerProduct(matrix_t::Ones(n, 1), B).eval();

    // Construct Pi by concatenating block matrices
    matrix_t Pi(C.cols() * Pi_1.rows(), Pi_1.cols() + Pi_2.cols());

    Pi << C(0, 0) * Pi_1, C(1, 0) * Pi_2,
          C(0, 1) * Pi_1, C(1, 1) * Pi_2,
          C(0, 2) * Pi_1, C(1, 2) * Pi_2,
          C(0, 3) * Pi_1, C(1, 3) * Pi_2;

    return Pi;
}

// Function to project onto the PSD cone
matrix_t Bezier::Proj_PSD(const matrix_t& M) {
    // ChatGPT
    Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(M);

    // Ensure the eigenvalue decomposition was successful
    if (eigensolver.info() != Eigen::Success) {
        throw std::runtime_error("Eigenvalue decomposition failed!");
    }

    // Get eigenvalues and eigenvectors
    vector_t evals = eigensolver.eigenvalues();
    matrix_t evecs = eigensolver.eigenvectors();

    // Initialize the output matrix with small positive values
    matrix_t M_ = matrix_t::Constant(M.rows(), M.cols(), std::numeric_limits<double>::epsilon());

    // Reconstruct the PSD matrix by adding contributions from positive eigenvalues
    for (int i = 0; i < evals.size(); ++i) {
        if (evals(i) > 0) {
            M_ += evals(i) * evecs.col(i) * evecs.col(i).transpose();
        }
    }

    return M_;
}

void Bezier::M_N_Gamma(const double Lg, const double Lf, const vector_t g_xbar, const double e_bar, const vector_t K, const double u_max,
                  matrix_t &M, vector_t &N, scalar_t &Gamma, vector_t &c, matrix_t &M_og) {
    // ChatGPT
    // Construct M_og
    M_og.resize(2, 2);
    M_og << Lg * Lf, Lg / 2.0,
            Lg / 2.0, 0.0;

    // Construct N
    N.resize(2);
    N << 2 * Lg * Lf * e_bar + Lf * g_xbar.norm() + Lg * K.norm() * e_bar,
         g_xbar.norm() + Lg * e_bar;

    // Compute Gamma
    Gamma = e_bar * (Lg * e_bar + g_xbar.norm()) * (Lf + K.norm());

    // Project M_og onto the PSD cone
    M = Proj_PSD(M_og);

    // Compute offset
    double offset = u_max - Gamma;

    // Solve quadratic equations for p1 and p2
    double p1 = (-N(0) + std::sqrt(N(0) * N(0) + 4 * M(0, 0) * offset)) / (2 * M(0, 0));
    double p2 = (-N(1) + std::sqrt(N(1) * N(1) + 4 * M(1, 1) * offset)) / (2 * M(1, 1));

    // Compute c based on eigenvalues
    c.resize(2);
    Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(M);
    if (eigensolver.eigenvalues().maxCoeff() < 1e-8) {
        c << 0, 1.0 / offset;
    } else {
        c << 1.0 / p1, 1.0 / p2;
    }
}

void Bezier::F_G(const matrix_t& Ax, const vector_t& bx, const matrix_t& H, const int m, const matrix_t& xbar, 
                const matrix_t& f_xbar, const matrix_t& g_xbar, const int gamma, const std::vector<matrix_t>& Q, 
                const double Lg, const double Lf, const double e_bar, const vector_t K, const double u_max,
                matrix_t &F, matrix_t &G) {
    // ChatGPT generated
    int order = H.rows();

    F.resize(0,H.cols());
    F.setZero();

    for (int j = 0; j < xbar.cols(); ++j) {
        // Compute H_vec and H_vec2
        matrix_t H_vec1 = H_vec(H, m, order - 1, gamma, gamma);
        matrix_t H_vec2 = H_vec(H, m, order - 1, gamma, gamma - 1);

        // Compute M, N, Gamma, and c using Bezier::M_N_Gamma
        matrix_t M, M_og;
        vector_t N, c;
        scalar_t Gamma;
        M_N_Gamma(Lg, Lf, g_xbar.block(0, j * m, g_xbar.rows(), m), e_bar, K, u_max,M,N,Gamma,c,M_og);

        // Compute Pi using Bezier::Pi
        matrix_t Pi_ = Pi(c, m, gamma);

        // Construct F
        matrix_t F_block = kroneckerProduct(matrix_t::Identity(order, order), Pi_) *
                           kroneckerProduct(Q[j].transpose(), kroneckerProduct(matrix_t::Identity(m, m), matrix_t::Identity(gamma + 1, gamma + 1))) *
                           H_vec1;
        F.conservativeResize(F.rows() + F_block.rows(), F.cols());
        F.block(F.rows() - F_block.rows(), F.cols() - F_block.cols(), F_block.rows(), F_block.cols()) = F_block;
        F_block = kroneckerProduct(matrix_t::Identity(order, order), Ax) *
                           kroneckerProduct(matrix_t::Identity(m, m), kroneckerProduct(Q[j].transpose(), matrix_t::Identity(gamma, gamma))) *
                           H_vec2;
        F.conservativeResize(F.rows() + F_block.rows(), F.cols());
        F.block(F.rows() - F_block.rows(), F.cols() - F_block.cols(), F_block.rows(), F_block.cols()) = F_block;

        // Construct G
        vector_t x_f(xbar.rows() + f_xbar.rows());
        x_f << xbar.col(j), f_xbar.col(j);
        vector_t G_block1 = kroneckerProduct(vector_t::Ones(order), vector_t::Ones(Pi_.rows()) + Pi_ * x_f);
        G.conservativeResize(G.rows() + G_block1.rows(), 1);
        G.block(G.rows() - G_block1.rows(), 0, G_block1.rows(), 1) = G_block1;

        vector_t G_block2 = kroneckerProduct(vector_t::Ones(order), bx);
        G.conservativeResize(G.rows() + G_block2.rows(), 1);
        G.block(G.rows() - G_block2.rows(), 0, G_block2.rows(), 1) = G_block2;
    }
}

matrix_t Bezier::K_matrix(int m, int n) {
    // ChatGPT generated
    int mn = m * n;
    matrix_t Y = matrix_t::Identity(mn, mn);

    // Generate the permutation indices
    Eigen::VectorXi I(mn);
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            I(i * n + j) = j * m + i;
        }
    }

    // Apply the permutation to the rows of the identity matrix
    matrix_t permuted_Y(mn, mn);
    for (int i = 0; i < mn; ++i) {
        permuted_Y.row(i) = Y.row(I(i));
    }

    return permuted_Y;
}


matrix_t Bezier::inv_DT_vec(int m, int order, int gamma) {
     matrix_t D = D_matrix(order, gamma);
     matrix_t kron_prod = kroneckerProduct(D.inverse().transpose(), matrix_t::Identity(m, m));
     return kron_prod;
}

matrix_t Bezier::R_n(int order1, int order2) {
  matrix_t R_big(order1+1, order2+1); 
  R_big.block(0,0,order1+1,order1+1).setIdentity();
  for (int r=order1+1; r<order2+1; r++) {
    matrix_t R_ = R_matrix(r);
    matrix_t R_int = R_big.block(0,0,order1+1,r)*R_.transpose();
    R_big.block(0,0,order1+1,r+1) << R_int;
  }
  return R_big;
}

matrix_t Bezier::D_matrix(int order, int gamma) {
  vector_t z1(order+1);
  vector_t z2(order+1);
  z1.setZero();
  z2.setZero();
  z1(0) = 1;
  z2(order) = 1;
  matrix_t D(order+1, 2*gamma);
  matrix_t H_(order+1,order+1);
  H_.setIdentity();
  for (int i = 0; i < gamma; i++) {
    D.block(0,i,order+1,1) << H_*z1;
    D.block(0,gamma+i,order+1,1) << H_*z2;
    H_ *= H;
  }
  return D;
}

matrix_t Bezier::M_matrix(int order) {
  matrix_t M(order+1, order+1);
  matrix_t P = Pascal_upper(order+1);
  matrix_t C(order+1, order+1);
  C.setZero();
  C.diagonal() << P.block(0,order,order+1,1);
  M = (P.partialPivLu().solve(C)).transpose();
  return M;
}

//https://rosettacode.org/wiki/Pascal_matrix_generation
matrix_t Bezier::Pascal_upper(int n) {
    matrix_t P(n,n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (j < i) P(i,j) = 0;
            else if (i == j || i == 0) P(i,j) = 1;
            else P(i,j) = P(i-1,j-1) + P(i,j-1);
        }
    }
    return P;
}


vector_t Bezier::T(scalar_t t) {
  vector_t T(order+1);
  for (int i = 0; i < order+1; i++) {
    T(i) = pow(t/tau,i);
  }
  return T; 
}

matrix_t Bezier::T(vector_t t) {
  matrix_t T(t.size(),order+1);
  for (int i = 0; i < order+1; i++) {
    T.block(0,i,t.size(),1) << (t/tau).array().pow(i);
  }
  return T; 
}

matrix_t Bezier::b(scalar_t t, matrix_t xi) {
  return T(t).transpose()*M*xi;
}

matrix_t Bezier::B(scalar_t t, matrix_t xi) {
  matrix_t B(1,gamma*xi.cols());
  for (int i = 0; i < gamma; i++) {
    B.block(0,i*xi.cols(),1,xi.cols()) << db(t, i, xi);
  }
  return B;
}

matrix_t Bezier::B(vector_t t, matrix_t xi) {
  matrix_t B(t.size(),gamma*xi.cols());
  for (int i = 0; i < gamma; i++) {
    B.block(0,i*xi.cols(),t.size(),xi.cols()) << db(t, i, xi);
  }
  return B;
}

matrix_t Bezier::b(vector_t t, matrix_t xi) {
  return T(t)*M*xi;
}

matrix_t Bezier::db(scalar_t t, int p, matrix_t xi) {
  Eigen::MatrixPower<matrix_t> Hpow(H);
  return T(t).transpose()*M*Hpow(p).transpose()*xi;
}

matrix_t Bezier::db(vector_t t, int p, matrix_t xi) {
  Eigen::MatrixPower<matrix_t> Hpow(H);
  return T(t)*M*Hpow(p).transpose()*xi;
}

// matrix_t Bezier::F(matrix_t A, vector_t b) {
//   
// }

void Bezier::mix_constraints(Eigen::Ref<matrix_t> A_mix, Eigen::Ref<vector_t> b_mix, vector_t c, vector_t x_bar, vector_t f_xbar) {
int n = 2;
int m = 1;

matrix_t A(n,n);
matrix_t B(m,n);

A.setZero();
A.block(m,0,n-m,n-m).setIdentity();
B.setZero();
B.block(0,n-m,m,m).setIdentity();

for (int i = 0; i < n; i++) {
  for (int j = 0; j < m; j++) {
    A_mix.block(4*i*m+j,0,4,n) << c(0)*A.block(i,0,1,n) + c(1)*B.block(j,0,1,n),
                            c(0)*A.block(i,0,1,n) - c(1)*B.block(j,0,1,n),
                            -c(0)*A.block(i,0,1,n) + c(1)*B.block(j,0,1,n),
                            -c(0)*A.block(i,0,1,n) - c(1)*B.block(j,0,1,n);
    b_mix.segment(4*i*m+j,4) << 1  +c(0)*x_bar(i) + c(1)*f_xbar(j),
                                1  +c(0)*x_bar(i) - c(1)*f_xbar(j),
                                1  -c(0)*x_bar(i) + c(1)*f_xbar(j),
                                1  -c(0)*x_bar(i) - c(1)*f_xbar(j);
  }
}
}

void Bezier::Pi_SD(matrix_t &A) {
   Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver(A);
   vector_t eval;
   matrix_t evec;
   eval = eigensolver.eigenvalues();
   evec = eigensolver.eigenvectors();
   A.setZero();
   for (int i = 0; i < eval.size(); i++) {
       if (eval(i) >= 0)
         A += eval(i)*evec.block(0,i,eval.size(),1)*evec.block(0,i,eval.size(),1).transpose();
   }
}

void Bezier::input_constraints(Eigen::Ref<matrix_t> A_u, Eigen::Ref<vector_t> b_u, vector_t c, vector_t x_bar, vector_t f_xbar) {
  Eigen::MatrixPower<matrix_t> Hpow(H);
  vector_t I_m(4);
  matrix_t A_tmp(2,4);
  matrix_t A(8,2);
  vector_t b(8);

  mix_constraints(A, b, c, x_bar, f_xbar);

  for (int i = 0; i < 4; i++) {
    I_m.setZero();
    I_m(i) = 1;
    A_tmp << I_m.transpose()*Hpow(0).transpose(), I_m.transpose()*Hpow(1).transpose();
    std::cout << A <<std::endl << std::endl;
    std::cout << A_tmp <<std::endl << std::endl;
    std::cout << D.transpose().inverse() <<std::endl << std::endl;
    std::cout << A*A_tmp*D.transpose().inverse() <<std::endl << std::endl;
    A_u.block(8*i,0,8,4) << A*A_tmp*D.transpose().inverse(); // TODO: optiimize
    b_u.segment(8*i, 8) << b;
  }

}



// void Bezier::vert2hyp(vector_t V, Eigen::Ref<matrix_t> A, Eigen::Ref<matrix_t> b) {
//  Eigen::Polyhedron poly; 
// }
// 
// matrix_t Bezier::hyp2vert(matrix_t A, matrix_t b) {
//   Eigen::Polyhedron poly;
//   bool success = poly.setHrep(A, b);
//   auto vrep = poly.vrep();
//   matrix_t V;
//   V = vrep.first;
//   return V;
// }
