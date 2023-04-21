#include <iostream>
#include <Eigen/Dense>
#include "../inc/Bezier.h"
#include "../inc/Types.h"
#include "Polyhedron.h"


int main() {
  Bezier B = Bezier(4,2,1);
  //std::cout << B.S << std::endl;
  //std::cout << B.R_matrix(4) << std::endl;
  //std::cout << B.R << std::endl;
  //std::cout << B.H << std::endl << std::endl;
  //
  //std::cout << B.R_n(4,7) << std::endl <<std::endl;

  //std::cout << B.D_matrix(4,2) <<std::endl << std::endl;

  //std::cout << B.M_matrix(5) << std::endl;
  //std::cout << B.T(0.5) << std::endl;
  vector_t xi(5);
  xi << 0,0,0,1,1;
  //std::cout << B.b(1, xi) << std::endl;
  matrix_t xi2(5,2);
  xi2 << 1,0,
         1,0,
         0,0,
         0,1,
         0,1;
  //std::cout <<xi2 <<std::endl<<std::endl;
  //std::cout << B.b(1, xi2) << std::endl<<std::endl;
  std::cout << B.B(0.5, xi2) << std::endl;
  vector_t t(5);
  t << 0,0.25, 0.5, 0.75, 1;
  //std::cout << B.b(t, xi) << std::endl;
  //std::cout << B.b(t, xi2) << std::endl;
  //std::cout << B.db(t,1, xi) << std::endl;

  Eigen::Polyhedron poly;
  matrix_t AHrepCone(4,2);
  vector_t bHrepCone(4);
  AHrepCone << 1, -1,
            -1, -1,
             -1, 1,
             1, 1;
  bHrepCone << 2, 2, 2, 2;
  bool success = poly.setHrep(AHrepCone, bHrepCone);
  auto vrep = poly.vrep();
  std::cout << std::endl<<std::endl;
  std::cout << vrep.first << std::endl;
}
