#include <Eigen/Dense>
#include "../inc/Types.h"

using namespace bezier;

class Bezier {
  public:
    int order;
    int gamma;
    scalar_t T;

    matrix_t S;
    matrix_t R;
    matrix_t H;
    matrix_t D;

    Bezier(int order, int gamma) {}

    Bezier(int order, int gamma, scalar_t T) : 
        order(order), gamma(gamma), T(T), 
        R(R_matrix(order)), S(S_matrix(order)), H(H_matrix(order)),
        D(D_matrix(order, gamma)){ }

    matrix_t H_matrix(int order);
    matrix_t S_matrix(int order);
    matrix_t R_matrix(int order);
    matrix_t R_n(int order1, int order2);
    matrix_t D_matrix(int order, int gamma);
};

int add(int i, int j) {
    return i + j;
}
