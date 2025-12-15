#include <iostream>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "gtest/gtest.h"



/*
  ____ The Matrix Class ____

  The three mandatory template parameters of Matrix are: 
  Matrix<typename Scalar, int RowsAtCompileTime, int ColsAtCompileTime>
  where:
  - Scalar is the scalar type, i.e. the type of the coefficients (int, float, double, ...)
  - RowsAtCompileTime and ColsAtCompileTime are the number of rows and columns of the matrix as known **at compile time**

  e.g.: th a N-by-M matrix of scalar type T, you would write 
  typedef Matrix<T, N, M> myNbyMmatrixofTypeT;

    typedef Eigen::Matrix<int,1,2> A;
    typedef Eigen::Matrix<int, 1, 2> B;

  Eigen is not limited to matrices whose dimensions are known at compile time. 
  The RowsAtCompileTime and ColsAtCompileTime template parameters can take the special value Dynamic which indicates that the size is unknown at compile time, so must be handled as a run-time variable.
  In Eigen terminology, such a size is referred to as a dynamic size; while a size that is known at compile time is called a fixed size.

  e.g.: to define a matrix of doubles with dynamic size, you would write
  typedef Matrix<double, Dynamic, Dynamic> MatrixXd;

*/


/** INSTRUCTIONS
 * follow the per-test instructions and fill-in their TODOs
 *
 * */


TEST(EigenTutorial, MatrixInitialization_4x4_integers){

  // Complete such that A4i is a 4-by-4 matrix of integers
  Eigen::Matrix<int,4,4/*TODO*/> A4i;
  ASSERT_EQ(A4i.rows(), 4);
  ASSERT_EQ(A4i.cols(), 4);
}

TEST(EigenTutorial, VectorInitialization_column_3_doubles){
  // Complete such that b3d is a column vector of 3 doubles?
  Eigen::Matrix<double, 3, 1/*TODO*/> b3d;
  ASSERT_EQ(b3d.rows(), 3);
  ASSERT_EQ(b3d.cols(), 1);
}

TEST(EigenTutorial, VectorInitialization_row_5_floats){
  //  a row vector of 5 floats?
  Eigen::Matrix<float, 1, 5/*TODO*/> bT5f;
  ASSERT_EQ(bT5f.rows(), 1);
  ASSERT_EQ(bT5f.cols(), 5);
}

TEST(EigenTutorial, MatrixInitialization_dynamic_3x3_doubles){
  // Complete to declare a 3x3 matrix of doubles
  Eigen::MatrixXd matrix3x3(3,3/*TODO*/);
  ASSERT_EQ(matrix3x3.rows(), 3);
  ASSERT_EQ(matrix3x3.cols(), 3);

}

TEST(EigenTutorial, VectorInitialization_OnTheFly_column_float){

  int N;
  std::cout << "Could you please set the number of rows for 'vectorNf' variable?" << std::endl;
  std::cin >> N;
  /*TODO*/Eigen::VectorXf vectorNf(N/*TODO*/);
  ASSERT_EQ(vectorNf.rows(), N);
  ASSERT_EQ(vectorNf.cols(), 1);
  
}
/*============================================================*/

TEST(EigenTutorial, MatrixAssignment_dense){

  Eigen::Matrix3i A3x3;
  Eigen::Matrix3i B3x3;
  // Fill A3x3 and B3x3 such that they are both
  // | 1 2 3 |
  // | 4 5 6 |
  // | 7 8 9 |
  // using two different methods:
  // 1. with the stream operator <<
  // 2. using per-element assignment operator (_,_)
  /* TODO */
  A3x3 << 1,2,3,
    4,5,6,
    7,8,9;

  int k(1);
  for(int i(0); i<3; i++){
      for(int j(0); j<3; j++){
          B3x3(i,j) = k;
          k++;
      }
  }
    /* end todo */

  ASSERT_EQ(A3x3.sum(), 45);
  ASSERT_EQ(A3x3.prod(), 362880);

  ASSERT_EQ(B3x3.sum(), 45);
  ASSERT_EQ(B3x3.prod(), 362880);

}


TEST(EigenTutorial, MatrixAssignment_sparse){

  using SintMat = Eigen::SparseMatrix<int>;// sparse matrix type
  using T = Eigen::Triplet<int>;// triplet, used to fill the sparse matrices
  
  SintMat sparseA;
  std::vector<T> triplets; // nonzero entries: T(i, j, a_ij) where the element at the i-th row and j-th column is nonzero and whose value is a_ij

  // Complete such that sparseA is a 5x5 sparse matrix equivalent to
  // | 0	3	0	0	0  |
  // | 22	0	0	0	17 |
  // | 7	5	0	1	0  |
  // | 0	0	0	0	0  |
  // | 0	0	14	0	8  |
  
  /* TODO */
  sparseA.resize(5,5);
  triplets.push_back(T(0,1,3));
  triplets.push_back(T(1,0,22));
  triplets.push_back(T(1,4,17));
  triplets.push_back(T(2,0,7));
  triplets.push_back(T(2,1,5));
  triplets.push_back(T(2,3,1));
  triplets.push_back(T(4,2,14));
  triplets.push_back(T(4,4,8));  
  /* end todo */

  sparseA.setFromTriplets(triplets.begin(), triplets.end());

  std::vector<std::vector<int>> fullA = {{0,3,0,0,0},
					 {22,0,0,0,17},
					 {7,5,0,1,0},
					 {0,0,0,0,0},
					 {0,0,14,0,8}};
  for (int k=0; k<sparseA.outerSize(); ++k)
    for (SintMat::InnerIterator it(sparseA,k); it; ++it) 
      ASSERT_EQ(it.value(), fullA[it.row()][it.col()]); 
  
}




int main(int _argc, char** _argv){

    testing::InitGoogleTest(&_argc, _argv);
    return RUN_ALL_TESTS();

}
