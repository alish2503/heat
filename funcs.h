#include <cmath> 

double sigma = 2;
double beta = sigma + 1;
double c = 5.0;
double delta = 1/8;

int n = 32;
double l = 1.0;
double T = 1.0;
double h, h2;
double coef = 2 * (sigma + 1) / (sigma*(sigma + 2));
double k0 = 0.01, q0 = 0.0001, te = 1.1;
double L_T = 2 * M_PI*sqrt(k0 / q0)*sqrt(sigma + 1) / sigma;
double eps = 1.e-3;
bool local_heat = true;

double u0(double t) {
    if(local_heat) return 0.0;
    return pow(2 * sigma * c * c * t, 1.0 / sigma);
}

double u1(double t) {
    if(local_heat) return 0.0;
    if(l <= c * t) return pow(2 * sigma * c * (c * t - l), 1.0 / sigma);
    return 0.0;
}

double f(double u) {
    if(local_heat) return q0 * pow(u, beta);
    return 0.0;
}
double k(double u) {
    if(local_heat) return k0 * pow(u, sigma);
    return 0.5 * pow(u, sigma);
}
double exact(double x, double t) {
    if(local_heat) {
        if(fabs(x) > L_T / 2) return 0.0;
        double cosin = cos((M_PI * x) / L_T)*cos((M_PI * x) / L_T);
        return pow(q0*(te - t), (-1.0 / sigma))*pow(coef*cosin, 1.0 / sigma);
    }
    if(x <= c * t) return pow(2 * sigma * c * (c * t - x), 1.0 / sigma);
    return 0.0;
}

