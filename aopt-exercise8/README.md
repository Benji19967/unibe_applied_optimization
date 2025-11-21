Applied Optimization course coding framework
=

This coding framework serves as a companion for the Applied Optimizaion course. 
You will have to fill functions throughout the successive exercises. 
In particular, each exercise will require you to fill specific parts of the code marked by the keyword "todo".


Introduction
-
First, you should know that this framework relies on CMake to build.
If it is not installed on your machine, you can find the relevant documentation here: https://cmake.org/install/

The main CMakeLists.txt is in the root folder and will build all the per-exercise projects, which are contained in dedicated sub-folder. 
New ones will be added over the weeks.

IDE
-
An Integrated Development Environment is basically a fancy code editor which allows you to do a lot of cool stuff and overall save you a lot of effort when coding projects.
You are **STRONGLY ENCOURAGED** to install such an IDE to work on this course's exercise.
There are several available depending on your platform and you are free to chose one.
However, we advise you to use JetBrain's CLion, which is free for students: https://www.jetbrains.com/clion/
Video guides on how to install and use CLion will come soon. 

Alternatives include (but are not limited to):
* QT Creator https://www.qt.io/product/development-tools
* Eclipse https://www.eclipse.org/ide/
* Visual Studio (Windows) https://visualstudio.microsoft.com/
* Xcode (Mac OS) https://developer.apple.com/xcode/


Matrices and Vectors types
-
The whole framework relies on Eigen (http://eigen.tuxfamily.org/index.php?title=Main_Page) to handle its mathematic types, namely Matrices and Vectors.
There are two style of Matrices:
* "Dense" Matrices, which are basically 2D vectors: https://eigen.tuxfamily.org/dox/group__TutorialMatrixClass.html
* "Sparse" Matrices:  which rely on an array of triplets, each triplet being (i, j, M_ij), i.e. contain the indices of a value inside the matrix and the value itself:  https://eigen.tuxfamily.org/dox/group__TutorialSparse.html
You should ** MOST DEFINITELY** start by reading both these pages to get familiar with those types.
Pay especially attention to the various constructors, initializers, accessors and standard operators provided by those classes. 

They are generally present in the course's code framework as "Mat", "Vec" and "SMat", type aliases for Eigen::MatrixXd, Eigen::VectorXd and Eigen::SparseMatrix<double> respectively.



FunctionBase, ParametricFunctionBase and FunctionBaseSparse
-
Those are the three basic classes representing functions. The idea is that all functions we'll consider in the coding exercises fall into either of the three categories.
That is why they all share the same interface consisting of:
* function evaluation
* gradient evaluation
* hessian matrix evaluation

The only difference between FunctionBase and FunctionBaseSpase is that the latter uses a Sparse Matrix to output its hessian.
The difference between Functionbase and ParametricFunctionBase is that the latter uses an additional parameter in its evaluation.

See the corresponding code files in include/FunctionBase/ for more details about the implementation



