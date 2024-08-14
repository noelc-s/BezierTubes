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

matrix_t Bezier::H_vec(matrix_t H, int m, int order, int gamma, int power_of_H) {
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
