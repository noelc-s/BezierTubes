#include <gtest/gtest.h>
#include "../inc/Bezier.h"
#include "../inc/Types.h"

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);

  int order = 4;
  int gamma = 2;
  scalar_t T = 1.;

  Bezier B = Bezier(order, gamma, T);
  matrix_t S_true(4,5);
  matrix_t R_true(5,4);
  matrix_t H_true(5,5);
  matrix_t R_n_true(5,8);
  matrix_t D_true(5,4);

  S_true << -4, 4, 0, 0, 0,
             0, -4, 4, 0, 0,
             0, 0, -4, 4, 0,
             0, 0, 0, -4, 4;
  R_true << 1, 0, 0, 0,
            0.25, 0.75, 0, 0,
            0, 0.5, 0.5, 0,
            0, 0, 0.75, 0.25,
            0, 0, 0, 1;
  H_true << -4, -1, 0, 0, 0,
            4, -2, -2, 0, 0,
            0, 3, 0, -3, 0,
            0, 0, 2, 2, -4,
            0, 0, 0, 1, 4;
  R_n_true <<         1,  0.428571,  0.142857, 0.0285714,         0,         0,         0,         0,
        0,  0.571429,  0.571429,  0.342857,  0.114286,         0,         0,         0,
        0,         0,  0.285714,  0.514286,  0.514286,  0.285714,         0,         0,
        0,         0,         0,  0.114286,  0.342857,  0.571429,  0.571429,         0,
        0,         0,         0,         0, 0.0285714,  0.142857,  0.428571,         1;
  D_true <<  1, -4,  0,  0,
        0,  4,  0,  0,
        0,  0,  0,  0,
        0,  0,  0, -4,
        0,  0,  1,  4;

  EXPECT_TRUE(B.S.isApprox(S_true,1e-6));
  EXPECT_TRUE(B.R.isApprox(R_true,1e-6));
  EXPECT_TRUE(B.H.isApprox(H_true,1e-6));
  EXPECT_TRUE(B.R_n(4,7).isApprox(R_n_true,1e-6));
  EXPECT_TRUE(B.D.isApprox(D_true,1e-6));
}

