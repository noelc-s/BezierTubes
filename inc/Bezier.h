#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "../inc/Types.h"
//#include "setoper.h"
//#include "cdd.h"
//#include "Polyhedron.h"

using namespace bezier;

class Bezier {
  public:
    int order;
    int gamma;
    scalar_t tau;

    matrix_t S;
    matrix_t R;
    matrix_t H;
    matrix_t D;
    matrix_t M;

    Bezier(int order, int gamma) {}

    Bezier(int order, int gamma, scalar_t tau) : 
        order(order), gamma(gamma), tau(tau), 
        R(R_matrix(order)), S(S_matrix(order)), H(H_matrix(order)),
        D(D_matrix(order, gamma)), M(M_matrix(order)){ }

    matrix_t H_matrix(int order);
    matrix_t M_matrix(int order);
    matrix_t S_matrix(int order);
    matrix_t R_matrix(int order);
    matrix_t R_n(int order1, int order2);
    matrix_t D_matrix(int order, int gamma);
    vector_t T(scalar_t tau);
    matrix_t T(vector_t tau);

    matrix_t b(scalar_t t, matrix_t xi);
    matrix_t b(vector_t t, matrix_t xi);
    matrix_t db(scalar_t t, int p, matrix_t xi);
    matrix_t db(vector_t t, int p, matrix_t xi);

    matrix_t B(scalar_t t, matrix_t xi);

  private:
    matrix_t Pascal_upper(int n);
};
