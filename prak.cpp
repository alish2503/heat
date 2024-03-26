#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <functional>
#include "funcs.h"

typedef std::vector<double> dv;

void non_linear_explicit_step(double, double, dv&, const dv&);
void implicit_solution(dv&, dv&, dv&, const dv&, const dv&, const dv&, const dv&, const dv&);
void non_linear_implicit_step(double, double, dv&, const dv&);
void iter_method(double, double, dv&, const dv&);
void plot();
void Runge_rule(double, int, dv&, const dv&, std::function<void(double, double, dv&, const dv&)>);
double error(const dv&, const dv&); 
double max(const dv&, const dv&);

int main(void) {
    double tau = 0.001;
    if(local_heat) l = 2 * L_T;
    h = l / n;
    h2 = h * h;
    dv y1(n + 1), y2(n + 1), x(n + 1), y_exact(n + 1);
    for (int i = 0; i <= n; i++) {
        x[i] = i * h;
        if(local_heat) x[i] -= L_T;
        y1[i] = y2[i] = phi2(x[i]);
    }
    std::ofstream ofs("points.txt");
    for (int i = 0; i <= n; i++) {
        y_exact[i] = exact(x[i], T);
        ofs << std::setprecision(15) << x[i] << ' ' << y_exact[i] << std::endl;
    }
    ofs << std::endl << std::endl;
    for (int i = 0; i <= n; i++) ofs << std::setprecision(30) << x[i] << ' ' << y1[i] << std::endl;
    ofs << std::endl << std::endl;
    ofs.close();
    Runge_rule(tau, 1, y1, x, &non_linear_explicit_step);
    Runge_rule(tau, 2, y2, x, &non_linear_implicit_step);
    if(local_heat) std::cout << L_T << std::endl;
    std::cout << "Explicit error: " << error(y1, y_exact) << std::endl;
    std::cout << "Implicit error: " << error(y2, y_exact) << std::endl;
    plot(); 
}

void Runge_rule(double tau, int num, dv &y, const dv &x, std::function<void(double, double, dv&, const dv&)> method) {
    dv y1, y2; y1 = y2 = y;
    double t = 0.0;
    char steps[] = "steps..txt"; steps[5] = char(num + 48);
    std::ofstream ofs1("points.txt", std::ios_base::app), ofs2(steps);
    while (true) {
        method(t, tau, y1, x); method(t, tau * 2, y2, x); 
        while (error(y1, y2) >= eps) {
            y2 = y1; y1 = y;
            tau /= 2;
            method(t, tau, y1, x); 
        }
        y = y2 = y1;
        t += tau;
        ofs2 << std::setprecision(6) << t << ' ' << tau << std::endl; 
        if(t >= T) break;
        y[0] = y1[0] = y2[0] = u0(t); 
        y[n] = y1[n] = y2[n] = u1(t);     
    }
    for (int i = 0; i <= n; i++) ofs1 << std::setprecision(6) << x[i] << ' ' << y[i] << std::endl; 
    ofs1 << std::endl << std::endl;
    ofs1.close(); ofs2.close();
}
void non_linear_explicit_step(double, double tau, dv &y, const dv &x) {
    double tmp = y[0];
    for (int i = 1; i < n; i++) {
        double yi = (0.5 * ((k(y[i]) + k(y[i + 1])) * (y[i + 1] - y[i]) - (k(y[i - 1]) + k(y[i])) * (y[i] - tmp)) / h2 + f(y[i])) * tau + y[i];
        tmp = y[i]; y[i] = yi;
    }
}
void non_linear_implicit_step(double t, double tau, dv &y, const dv &x) {
    dv alpha(n), beta(n), A(n), B(n), C(n), F(n);
    double gamma = tau / h2;
    double next_tj = t + tau;
    for (int i = 1; i < n; i++) { 
        A[i] = gamma * 0.5 * (k(y[i - 1]) + k(y[i]));
        C[i] = gamma * 0.5 * (k(y[i]) + k(y[i + 1]));
        B[i] = -(1.0 + A[i] + C[i]);
    }
    for (int i = 2; i < n - 1; i++) F[i] = -(tau * f(y[i]) + y[i]);
    F[1] = -(tau * f(y[1]) + y[1]) - A[1] * u0(next_tj);
    F[n - 1] = -(tau * f(y[n - 1]) + y[n - 1]) - C[n - 1] * u1(next_tj);
    implicit_solution(alpha, beta, y, x, A, B, C, F);
}
void implicit_solution(dv &alpha, dv &beta, dv &y, const dv &x, const dv &A, const dv &B, const dv &C, const dv &F) {
    alpha[2] = -C[1] / B[1];
    beta[2] = F[1] / B[1];
    for (int i = 2; i < n - 1; i++) {  
        double tmp = (A[i] * alpha[i] + B[i]);
        alpha[i + 1] = -C[i] / tmp;
        beta[i + 1] = (F[i] - A[i] * beta[i]) / tmp;
    }
    y[n - 1] = (F[n - 1] - A[n - 1] * beta[n - 1]) / (B[n - 1] + A[n - 1] * alpha[n - 1]);
    for (int i = n - 2; i > 0; i--) y[i] = alpha[i + 1] * y[i + 1] + beta[i + 1];
}
void iter_method(double t, double tau, dv &y, const dv &x) {
    dv y_next = y;
    non_linear_implicit_step(t, tau, y_next, x);
    int iter = 0;
    while (max(y, y_next) >= eps) { 
        y = y_next;
        if(iter == 2) break;
	non_linear_implicit_step(t, tau, y_next, x);
        iter++;
	}
    y = y_next; 
}
double error(const dv &y1, const dv &y2) {
	double error = 0.0;
	for (int i = 1; i < n; i++) error += (y1[i] - y2[i])*(y1[i] - y2[i]);
    error *= h;	
	return sqrt(error);
}
double max(const dv &y1, const dv &y2) {
	double max = fabs(y1[1] - y2[1]);
	for (int i = 2; i < n; i++) {
        double dif = fabs(y1[i] - y2[i]);
		if(max < dif) max = dif;
    }
	return max;
}
void plot() {
    FILE *gp = popen("gnuplot","w"); 
	if(!gp) {
        printf("Error opening pipe to GNU plot\n");
        return;
    }
    if(local_heat) fprintf(gp, "set output 'local_heat_time_steps.png'\n"); 
    else fprintf(gp, "set output 'problem_time_steps.png'\n");
    fprintf(gp, "set terminal png\n");
    fprintf(gp, "plot 'steps1.txt' lt rgb 'blue' title 'explicit','steps2.txt' lt rgb 'green' title 'implicit'\n");   
    fprintf(gp, "set xlabel 'x'\n");
    fprintf(gp, "set ylabel 'y'\n");
    fprintf(gp, "set terminal png\n");
    if(local_heat) {
        fprintf(gp, "set output 'local_heat.png'\n");
        fprintf(gp, "set label 'length: %lf' at 30,230\n", L_T);
    }
    else fprintf(gp, "set output 'problem.png'\n");
    fprintf(gp, "plot 'points.txt' index 0 with linespoints pt 7 lt rgb 'green' title 'exact','points.txt' index 1 with linespoints pt 7 lt rgb 'blue' title 'init','points.txt' index 2 with linespoints pt 7 lt rgb 'orange' title 'explicit','points.txt' index 3 with linespoints pt 7 lt rgb 'yellow' title 'implicit'\n"); 
    fclose(gp);
    remove("points.txt"); remove("steps1.txt"); remove("steps2.txt");
}

