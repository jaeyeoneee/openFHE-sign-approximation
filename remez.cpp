#include <iostream>
#include <vector>

void ChebyshevBasis(int n, double x, double w, std::vector<double>& rs) {
    rs[0] = 1;
    if (n == 0) {
        return;
    } 

    rs[1] = x/w;
    if (n == 1) {
        return;
    }

    auto compute_T = [&](int i, int j) {
        rs[i+j] = 2 * rs[j] * rs[i] - rs[j-i];
    };

    for (int k = 2; k <= n; k++) {
        int i = k / 2;
        int j = (k + 1) / 2;
        compute_T(i, j);
    }
}

auto ImprovedRemez(int n, int domain[], double gamma, int iter, double w){
    // TODO variable 수정
    // Step 1: Choose initial points
    

} 

int main() {
    // Chebyshev polynomials test
    int cn = 5;
    double cx = 0.5, cw = 2;
    std::vector<double> rs(cn + 1);
    ChebyshevBasis(cn, cx, cw, rs);

    std::cout << "<<------Chebyshev polynomials------>>" << std::endl;
    for (int i = 0; i <= cn; i++) {
        std::cout << "rs" << i << " = " << rs[i] << std::endl;
    }

    //Algorithm2: Improved Multi-Interval Remes Algorithm test
    int domain[2] = {0.1, 1};
    int n = 5, iter = 10;
    double gamma = 0.1, w = 2;
    ImprovedRemez(n, domain, gamma, iter, w);

    return 0;
}