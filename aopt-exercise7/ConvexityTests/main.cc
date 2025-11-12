#include <iostream>
#include <Algorithms/ConvexityTest.hh>
#include <Functions/FunctionNonConvex2D.hh>
#include <Functions/FunctionQuadraticND.hh>


void testNonConvexFunction() {
    // try different domains if you want
    AOPT::FunctionNonConvex2D function;
    std::cout << "Checking the convexity of a 2-D function ..." << std::endl;
    AOPT::ConvexityTest::isConvex(&function);
}

void testQuadraticND(int n) {
    AOPT::FunctionQuadraticND function(n);
    std::cout << "Checking the convexity of a "<<n<<"-D function ..." << std::endl;
    AOPT::ConvexityTest::isConvex(&function);
}

int main(int argc, const char* argv[] ) {
    if (argc < 2) {
        std::cout << "usage:\n" << argv[0] << " option [dimension=3]\n"
                  << "1. test Nonconvex 2d function\n"
                  << "2. test QuadraticND function\n";

        return -1;
    }
    int opt = atoi(argv[1]);
    int n = 3;
    if (argc > 2) n = atoi(argv[2]);
    switch (opt) {
        case 1:
            testNonConvexFunction();
            break;
        case 2:
            testQuadraticND(n);
            break;
        default:
            std::cout << "option " << opt << " not handled by any case" << std::endl;
    }
    return 0;
}
