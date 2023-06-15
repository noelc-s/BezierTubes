#include <iostream>
#include <Eigen/Dense>
#include "../inc/Bezier.h"
#include "../inc/Types.h"
#include "Polyhedron.h"


int main() {
  // Arguments of Bezier are: order of curve (must be at least 3), relative degree of output (2), 
  // and time over which the curve is to be defined
  Bezier B = Bezier(3,2,1);
  // These are the control points. In this example they are a vector:
  vector_t xi(4); // of size (order+1)
  xi << 0,0,1,1;
  // This evaluates the bezier curve at the time (the first argument) 1.
  std::cout << B.b(1, xi) << std::endl;
  // The control points can also be a matrix (for multiple outputs, e.g. x,y
  matrix_t xi2(4,2);
  xi2 << 1,0,
         1,0,
         0,1,
         0,1;
  std::cout << B.b(1, xi2) << std::endl<<std::endl;
  vector_t t(5);
  t << 0,0.25, 0.5, 0.75, 1;
  // The above functions can also be used with a sequence of time points
  std::cout << B.b(t, xi) << std::endl;
  std::cout << B.b(t, xi2) << std::endl;
  // This is the bezier curve in state space evaluated at the time points t
  std::cout << B.B(t, xi2) << std::endl;
  // And this is the time derivative
  std::cout << B.db(t,1, xi) << std::endl;
  //
  ///////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  ///std::cout << B.S << std::endl;
  //std::cout << B.R_matrix(4) << std::endl;
  //std::cout << B.R << std::endl;
  //std::cout << B.H << std::endl << std::endl;
  //
  //std::cout << B.R_n(4,7) << std::endl <<std::endl;

  //std::cout << B.D_matrix(4,2) <<std::endl << std::endl;

  //std::cout << B.M_matrix(5) << std::endl;
  //std::cout << B.T(0.5) << std::endl;
/

  scalar_t Lf, Lg;
  vector_t x_bar(2);
  vector_t f_xbar(1);
  vector_t g_xbar(1);
  scalar_t offset;
  x_bar << 0,0;
  f_xbar << 0;
  g_xbar << .1;
  offset = 1;
  Lf = 1;
  Lg = 1;

  matrix_t M(2,2);
  vector_t N(2);
  M << Lg*Lf, Lg/2, Lg/2, 0;
  B.Pi_SD(M);
  //std::cout << M << std::endl;
  N << Lf*g_xbar.norm(), g_xbar.norm();
  //std::cout << N << std::endl;

  scalar_t p1, p2;
  p1 = (-N(0) + sqrt(pow(N(0),2)+4*M(0,0)*offset))/(2*M(0,0));
  p2 = (-N(1) + sqrt(pow(N(1),2)+4*M(1,1)*offset))/(2*M(1,1)); 

  Eigen::Polyhedron poly;

  int n = 2;
  int m = 1;
  vector_t c(2);
  c << 1./p1,1./p2;
  std::cout << c << std::endl;
  matrix_t A_mix(32,4);        
  vector_t b_mix(32);
  B.input_constraints(A_mix, b_mix, c, x_bar, f_xbar);
  std::cout << A_mix << std::endl<<std::endl;
  std::cout << b_mix << std::endl<<std::endl;
  matrix_t V;
  V = B.hyp2vert(A_mix.block(0,0,32,2), b_mix);
  //bool success = poly.setHrep(A_mix, b_mix);
  //auto vrep = poly.vrep();
  //std::cout << std::endl<<std::endl;
  //std::cout << vrep.first << std::endl;
  //std::cout << vrep.second << std::endl;
  std::cout <<  V << std::endl;
}
