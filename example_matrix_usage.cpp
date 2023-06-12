#include <iostream>
#include "matrix.h"

int main() {
    /// @brief this code demonstrates using Vec and Mat classes

    
    /// @brief declaration of vectors
    std::cout << "---Declaration of vectors---\n";
    Vec<double, 2> p0;              // 2-dim vector, each component has `double` type, initialize to 0
    Vec<double, 3> p1(1, 2, 4.2);   // 3-dim vector, each component has `double` type, initialize to given values
    Vec<int, 4> p2(2, 2, 7, 4);     // 4-dim vector, each component has `int` type, initialize to given values

    // you can print them with cout
    std::cout << "p0 = " << p0 << "\n";
    std::cout << "p1 = " << p1 << "\n";
    std::cout << "p2 = " << p2 << "\n";
    
    std::cout << "\n---Access components of vectors---\n";
    std::cout << "p1[2] = " << p1[2] << "\n";
    p1[2] = 0.9;
    std::cout << "setting p1[2] to 0.9, and now: p1 = " << p1 << "\n";



    /// @brief vector arithmetrics
    std::cout << "\n---Arithmetics of vectors---\n";
    Vec<double, 3> v1(2, 0, 1);
    Vec<double, 3> v2(1, 3, 1);
    double a = 2;

    std::cout << "v1 = " << v1 << "\n";
    std::cout << "v2 = " << v2 << "\n";
    std::cout << "a  = " << a  << "\n";
    std::cout << "v1 + v2 = " << v1 + v2 << "\n";       // vector addition
    std::cout << "v1 - v2 = " << v1 - v2 << "\n";       // vector subtraction
    std::cout << "v1 * v2 = " << v1 * v2 << "\n";       // elementwise multiplication
    std::cout << "v1 * a  = " << v1 * a  << "\n";       // scalar multiplication
    std::cout << "v1 / a  = " << v1 / a  << "\n";       // scalar multiplication
    std::cout << "<v1,v2> = " << v1.dot(v2) << "\n";    // dot product
    std::cout << "v1 x v2 = " << v1.cross(v2) << "\n";  // cross product

    std::cout << "square norm(v1) = " << v1.square_norm() << "\n";   // square norm
    std::cout << "v1 / norm(v1) = " << v1.normalized() << "\n";   // normalization

    /// @brief declaration of matrices
    std::cout << "\n---Declaration of matrices---\n";
    Mat<double, 3, 4> m1;   // a 3-by-4 matrix, initialized to 0 by default
    Mat<double, 3, 4> m2{3, 7, 7, 5, 3, 2}; // initialize with an initializer list (row-major indexing), if there are less than 12 numbers, the rest are filled with 0
    Mat<double, 3, 4> m3{2, -1, 0, 2, 4, 2}; // initialize with an initializer list (row-major indexing), if there are less than 12 numbers, the rest are filled with 0
    Mat<double, 4, 2> m4{1, -1, 0, 2, 7, 7, 5, 3, 2}; // if there are more than 8 numbers, the extra ones are ignored.
    Vec<double, 4> w(0, 1, 0, 1);
    std::cout << "m1 =\n" << m1 << "\n";
    std::cout << "m2 =\n" << m2 << "\n";
    std::cout << "m3 =\n" << m3 << "\n";
    std::cout << "m4 =\n" << m4 << "\n";
    std::cout << "w  = " << w << "\n";

    std::cout << "\n---Access items of matrices---\n";
    std::cout << "m2(2,1) = " << m2(2,1) << "\n";
    m2(2, 1) = 0.9;
    std::cout << "setting m2(2,1) to 0.9, and now: m2 =\n" << m2 << "\n";


    /// @brief matrix arithmetrics
    std::cout << "\n---Arithmetics of matrices---\n";
    std::cout << "m2 + m3 =\n" << m2 + m3 << "\n\n";  // addition
    std::cout << "m2 - m3 =\n" << m2 - m3 << "\n\n";  // subtraction
    std::cout << "m2 * m3 =\n" << m2 * m3 << "\n\n";  // elementwise multiplication
    std::cout << "matmul(m3, m4) =\n" << m3.matmul(m4) << "\n\n";  // matrix multiplication
    std::cout << "matmul(m3, w) = "  << m3.matmul(w)  << "\n\n";   // matrix-vector multiplication
    std::cout << "m3.transpose() =\n" << m3.transpose() << "\n\n";  // addition

    /// @note determinants (supports only matrices smaller than 3x3)
    Mat<double, 3, 3> m5{2, -1, 0, 2, 4, 2, 1, 1, 1}; // initialize with an initializer list (row-major indexing), if there are less than 12 numbers, the rest are filled with 0
    std::cout << "m5 =\n" << m5 << "\n";
    std::cout << "det(m5) = " << m5.det() << "\n";  // addition
}