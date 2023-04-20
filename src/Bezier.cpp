#include "../inc/Bezier.h"

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
  H_matrix = (1./T)*S.transpose()*R.transpose();
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
