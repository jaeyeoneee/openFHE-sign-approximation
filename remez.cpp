#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>


double f(double x){
    if (x == 0) return 0;
    else if (x > 0) return 1;
    else return -1; 
}


std::vector<double> ChebyshevBasis(int n, double x, double w) {
    std::vector<double> rs(n);
    n--;

    rs[0] = 1;
    if (n == 0) {
        return rs;
    } 

    rs[1] = x/w;
    if (n == 1) {
        return rs;
    }

    auto compute_T = [&](int i, int j) {
        rs[i+j] = 2 * rs[j] * rs[i] - rs[j-i];
    };

    for (int k = 2; k <= n; k++) {
        int i = k / 2;
        int j = (k + 1) / 2;
        compute_T(i, j);
    }

    return rs;
}


std::vector<double> SelectInitialNode(int n, double domain[]){
    double a = domain[0], b = domain[1];
    std::vector<double> nodes(n+1);

    for (int i=1; i<=n+1; i++){
        nodes[i-1] = (a+b)/2 + (b-a)/2 * std::cos((2.0*i-1)/(2.0*n+2)*M_PI);
    }

    return nodes;
}


std::vector<double> ComputeCoefficents(int n, std::vector<double> nodes, double w){
    Eigen::MatrixXd A(n+1, n+1);
    Eigen::VectorXd b(n+1);

    for (int i=0; i<n+1; i++){
        std::vector<double> basis = ChebyshevBasis(n, nodes[i], w);
        for (int k=0; k<n; k++){
            A(i, k) = basis[k];
        }
        A(i, n) = ((i+1)%2 == 0) ? 1 : -1;
        b(i) = f(nodes[i]);
    }

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    std::vector<double> coeffs(n+1);
    for (int i = 0; i < n+1; i++) {
        coeffs[i] = x(i);
    }

    return coeffs;
}

                                                                                                                                                                  
void ImprovedRemez(int n, double domain[], double gamma, int iter, double w){
    // TODO variable 수정
    // Step 1: Choose initial points
    std::vector<double> rs = SelectInitialNode(n, domain);

    // Step 2: find the polynomial p(x) in terms of chebyshev basis

} 

int main() {
    // variables
    int n = 5, iter = 10;
    double gamma = 0.1, w = 2;
    double domain[2] = {0.1, 1};
    
    // Chebyshev polynomials test
    double cx = 0.5;
    std::vector<double> rs = ChebyshevBasis(n, cx, w);
    std::cout << "<<------Chebyshev polynomials------>>" << std::endl;
    std::cout << "Chebyshev polynomials at x = " << cx << std::endl;
    for (int i = 0; i < n; i++) {
        std::cout << "rs" << i << " = " << rs[i] << std::endl;
    }
    std::cout << std::endl;

    //SelectInitialNode test
    std::vector<double> nodes = SelectInitialNode(n, domain);
    std::cout << "<<------SelectInitialNode------>>" << std::endl;
    for (int i = 0; i <= n; i++) {
        std::cout << "node" << i << " = " << nodes[i] << std::endl;
    }
    std::cout << std::endl;

    //ComputeCoefficents test
    std::vector<double> coeffs = ComputeCoefficents(n, nodes, w);
    std::cout << "<<------ComputeCoefficents------>>" << std::endl;
    for (int i = 0; i <= n; i++) {
        std::cout << "coeff" << i << " = " << coeffs[i] << std::endl;
    }

    //Algorithm2: Improved Multi-Interval Remes Algorithm test
    ImprovedRemez(n, domain, gamma, iter, w);

    return 0;
}