#include <iostream>
#include <Eigen/Dense>
#include "../inc/Bezier.h"
#include "../inc/Types.h"


int main() {
  Bezier B = Bezier(4,2,1);
  std::cout << B.S << std::endl;
  std::cout << B.R_matrix(4) << std::endl;
  std::cout << B.R << std::endl;
  std::cout << B.H << std::endl << std::endl;
  
  std::cout << B.R_n(4,7) << std::endl <<std::endl;

  std::cout << B.D_matrix(4,2) <<std::endl;


}
