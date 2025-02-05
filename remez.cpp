#include "Point.h"
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>

// Target function (sign function)
double f(double x){
    if (x == 0) return 0;
    else if (x > 0) return 1;
    else return -1; 
}


// Step 1: Select Initial Chebyshev Nodes
std::vector<double> SelectInitialNode(int n, double domain[]){
    double a = domain[0], b = domain[1];
    std::vector<double> nodes(n+2);

    for (int i=0; i<n+2; i++){
        nodes[i] = (a+b)/2 + (b-a)/2 * std::cos((2.0*(i+1)-1)/(2.0*(n+1)+2)*M_PI);
    }

    return nodes;
}


// Compute Chebyshev Basis
std::vector<double> ChebyshevBasis(int n, double x, double w) {
    std::vector<double> rs(n+1);

    rs[0] = 1;
    if (n == 0) {
        return rs;
    } 

    rs[1] = x/w;
    if (n == 1) {
        return rs;
    }

    for (int k = 2; k<n+1; k++){
        rs[k] = 2*x/w*rs[k-1] - rs[k-2];
    }

    return rs;
}


// Step 2:  Compute Polynomial Coefficients and Maximum Error
std::vector<double> ComputeCoefficents(int n, std::vector<double> nodes, double w){
    Eigen::MatrixXd A(n+2, n+2);
    Eigen::VectorXd b(n+2);

    for (int i=0; i<n+2; i++){
        std::vector<double> basis = ChebyshevBasis(n, nodes[i], w);
        for (int k=0; k<n+1; k++){
            A(i, k) = basis[k];
        }
        A(i, n+1) = ((i+1)%2 == 0) ? 1 : -1;
        b(i) = f(nodes[i]);
    }

    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    std::vector<double> coeffs(n+2);
    for (int i = 0; i < n+2; i++) {
        coeffs[i] = x(i);
    }

    return coeffs;
}

// Evaluate Chebyshev Polynomial
double ChebyshevEvalError(int n, const std::vector<double>& coeffs, double x, double w){
    std::vector<double> rs = ChebyshevBasis(n, x, w); 
    double result = 0;
    for (int i = 0; i < n+1; i++){
        result += coeffs[i]*rs[i];
    }
    result -= f(x);
    return result;
}

// Step 3: Collect Extreme Points 
std::vector<Point> CollectExtremePoints(int n, double sc, double precision, double domain[], const std::vector<double>& coeffs, double w){
    double a = domain[0], b = domain[1];
    std::vector<Point> extreme_points;

    double xPrev1 = a;
    double xPrev2 = a + sc;
    double xCurr = a + 2*sc;
    double yPrev1 = ChebyshevEvalError(n, coeffs, xPrev1, w);
    double yPrev2 = ChebyshevEvalError(n, coeffs, xPrev2, w);

    while (xCurr < b){
        double yCurr = ChebyshevEvalError(n, coeffs, xCurr, w);


        if ((yPrev1>yPrev2 && yPrev1>yCurr)||(yPrev1<yPrev2 && yPrev1<yCurr)){
            double left = xPrev2;
            double right = xCurr;
            bool isMax = yPrev1 > yPrev2;

            for (int i=0; i<precision; i++){
                double mid = (left+right)/2;

                double delta = (right-left)/4;
                std::vector<double> xPoints = {mid - delta, mid, mid + delta};
                std::vector<double> yPoints = {
                    ChebyshevEvalError(n, coeffs, xPoints[0], w), 
                    ChebyshevEvalError(n, coeffs, xPoints[1], w), 
                    ChebyshevEvalError(n, coeffs, xPoints[2], w)
                };

                int argmax;
                if (isMax){
                    argmax = std::distance(yPoints.begin(), std::max_element(yPoints.begin(), yPoints.end()));
                }
                else {
                    argmax = std::distance(yPoints.begin(), std::min_element(yPoints.begin(), yPoints.end()));
                }

                if (argmax == 0) right = mid;
                else if (argmax == 2) left = mid;
                else {
                    left = xPoints[0];
                    right = xPoints[2];
                }
            }

            double xExtreme = (left+right)/2;
            double yExtreme = ChebyshevEvalError(n, coeffs, xExtreme, w);
            
            if (isMax && yExtreme > std::abs(coeffs[n+1])){
                extreme_points.push_back(Point(xExtreme, yExtreme, (isMax) ? 1 : -1));
            }
            else if (!isMax && (-1)*yExtreme > std::abs(coeffs[n+1])){
                extreme_points.push_back(Point(xExtreme, yExtreme, (isMax) ? 1 : -1));
            }
        }

        xPrev2 = xPrev1;
        xPrev1 = xCurr;
        yPrev2 = yPrev1;
        yPrev1 = yCurr;
        xCurr += sc;
    }

    return extreme_points;
}

class RefPoint {
    public:
        double add;
        int ind;

        RefPoint(double _add, int _ind){
            add = _add;
            ind = _ind;
        }
};

// Step 4: Find n+1 extreme points
std::vector<Point> FindExtremePoints(int n, double domain[], std::vector<Point>& extremes ,const std::vector<double>& coeffs, double w){

    if (extremes.size() < n+ 2){
        throw std::runtime_error("Not enough extreme points for the given degree");
    }

    std::cout << "extremes.size() before removing points = " << extremes.size() << std::endl;

    // Step 3-1~8: Remove points that don't satisfy alternating condition
    int i = 0;
    while (i < extremes.size() - 1){
        if (extremes[i].locmm * extremes[i+1].locmm == 1){
            if (std::abs(extremes[i].y) >= std::abs(extremes[i+1].y)){
                extremes.erase(extremes.begin() + i+1);
            }
            else {
                extremes.erase(extremes.begin() + i);
            }
        }
        else {
            i++;
        }
    }
    
    std::cout << "extremes.size() after alternating condition = " << extremes.size() << std::endl;

    // Step 3-9~24: If more than d+2 points remain, remove points to get exactly d+2
    // TODO:순서 변경
    while (extremes.size() > n+2){
        if (extremes.size() == n + 3){
            if (std::abs(extremes[0].y) >= std::abs(extremes[extremes.size()-1].y)){
                extremes.pop_back();
            }
            else {
                extremes.erase(extremes.begin());
            }
        }
        else {
            std::vector<RefPoint> B;
            for (int i = 0; i < extremes.size()-1; i++){
                B.push_back(RefPoint(std::abs(extremes[i].y) + std::abs(extremes[i+1].y), i));
            }
            if (extremes.size() == n+4){
                B.push_back(RefPoint(std::abs(extremes[extremes.size()-1].y)+std::abs(extremes[0].y), -1));
                std::sort(B.begin(), B.end(), [](RefPoint a, RefPoint b) {return a.add < b.add;});
                int ind = B[0].ind;
                if (ind == -1){
                    extremes.erase(extremes.begin());
                    extremes.pop_back();
                }
                else {
                    extremes.erase(extremes.begin() + ind);
                    extremes.erase(extremes.begin() + ind+1);
                }
            }
            else {
                std::sort(B.begin(), B.end(), [](RefPoint a, RefPoint b) {return a.add < b.add;});
                int ind = B[0].ind;
                if (ind == 0){
                    extremes.erase(extremes.begin());
                }
                else if (ind == extremes.size() - 2){
                    extremes.pop_back();
                }
                else {
                    extremes.erase(extremes.begin() + ind);
                    extremes.erase(extremes.begin() + ind+1);
                }
            }

        }
    }

    std::cout << "extremes.size() after removing points = " << extremes.size() << std::endl;

    return extremes;
}
                                                                                                                                                                  
std::vector<double> ImprovedRemez(int n, double domain[], double gamma, int iter, double w, double sc, double precision, double epsilon){
    // Step 1: Choose initial points
    std::vector<double> rs = SelectInitialNode(n, domain);

    while (true){
        // Step 2: find the polynomial p(x) in terms of chebyshev basis
        std::vector<double> coeffs = ComputeCoefficents(n, rs, w);

        // Step 3: Collect extreme points
        std::vector<Point> extreme_points = CollectExtremePoints(n, sc, precision, domain, coeffs, w);

        // Step 4: Find n+2 extreme points
        std::vector<Point> extremes = FindExtremePoints(n, domain, extreme_points, coeffs, w);

        auto max_iter = std::max_element(extremes.begin(), extremes.end(), 
        [](const Point& a, const Point& b) { return std::abs(a.y) < std::abs(b.y); });

        auto min_iter = std::min_element(extremes.begin(), extremes.end(), 
        [](const Point& a, const Point& b) { return std::abs(a.y) < std::abs(b.y); });

        double max_error = (max_iter != extremes.end()) ? std::abs(max_iter->y) : 0;
        double min_error = (min_iter != extremes.end()) ? std::abs(min_iter->y) : 0;

        std::cout << "max_error = " << max_error << std::endl;
        std::cout << "min_error = " << min_error << std::endl;

        if ((max_error - min_error)/min_error < epsilon){
            return coeffs;
        }

        for (int i=0; i<n+2; i++){
            rs[i] = extremes[i].x;
        }
    }
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
    for (int i = 0; i < n+1; i++) {
        std::cout << "rs" << i << " = " << rs[i] << std::endl;
    }
    std::cout << std::endl;

    //SelectInitialNode test
    std::vector<double> nodes = SelectInitialNode(n, domain);
    std::cout << "<<------SelectInitialNode------>>" << std::endl;
    for (int i = 0; i < n+2; i++) {
        std::cout << "node" << i << " = " << nodes[i] << std::endl;
    }
    std::cout << std::endl;

    //ComputeCoefficents test
    std::vector<double> coeffs = ComputeCoefficents(n, nodes, w);
    std::cout << "<<------ComputeCoefficents------>>" << std::endl;
    for (int i = 0; i < n+1; i++) {
        std::cout << "coeff" << i << " = " << coeffs[i] << std::endl;
    }
    std::cout << "error" << " = " << coeffs[n+1] << std::endl;

    //ChebyshevEvalError test
    double result = ChebyshevEvalError(n, coeffs, cx, w);
    std::cout << "<<------ChebyshevEvalError------>>" << std::endl;
    std::cout << "ChebyshevEvalError at x = " << cx << " = " << result << std::endl;

    //CollectExtremePoints test
    double sc = 0.000001, precision = 10;
    std::vector<Point> extreme_points = CollectExtremePoints(n, sc, precision, domain, coeffs, w);
    std::cout << "<<------CollectExtremePoints------>>" << std::endl;
    for (int i = 0; i < extreme_points.size(); i++) {
        std::cout << "extreme_point" << i << " = " << extreme_points[i].x << ", " << extreme_points[i].y << "," << extreme_points[i].locmm << std::endl;
    }

    //FindExtremePoints test
    std::cout << std::endl;
    std::cout << "<<------FindExtremePoints------>>" << std::endl;
    std::vector<Point> extremes = FindExtremePoints(n, domain, extreme_points, coeffs, w);
    for (int i = 0; i < extremes.size(); i++) {
        std::cout << "extreme_point" << i << " = " << extremes[i].x << ", " << extremes[i].y << "," << extremes[i].locmm << std::endl;
    }    

    //Algorithm2: Improved Multi-Interval Remes Algorithm test
    std::cout << "<<------Improved Multi-Interval Remes Algorithm------>>" << std::endl;
    double epsilon = 0.0000001;
    std::vector<double> coeffRs = ImprovedRemez(n, domain, gamma, iter, w, sc, precision, epsilon);

    return 0;
}