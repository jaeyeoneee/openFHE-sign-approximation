#include <iostream>
#include <vector>

void chebyshev (int n, double x, double w, std::vector<double>& rs) {
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



int main() {
    // Chebyshev polynomials test
    int n = 5;
    double x = 0.5, w = 2;
    std::vector<double> rs(n + 1);
    chebyshev(n, x, w, rs);

    std::cout << "<<------Chebyshev polynomials------>>" << std::endl;
    for (int i = 0; i <= n; i++) {
        std::cout << "rs" << i << " = " << rs[i] << std::endl;
    }

    //Algorithm2: Improved Multi-Interval Remes Algorithm test



    return 0;
}