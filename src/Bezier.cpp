#include "../inc/Bezier.h"
#include <iostream>

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
