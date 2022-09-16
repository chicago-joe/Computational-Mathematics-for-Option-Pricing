//  Computational Homework #3
//  UIUC - IE525 - Spring 2019
//
//  Created by Joseph Loss on 4/10/2019
//
// Question 3

#include "distributions.hpp"
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <chrono>
#include <random>
#include <vector>
#include <fstream>
#include <string>
using namespace std;

#define E 2.718281828459045

#ifndef  PI
const double PI = 3.141592653589793238462643;
#endif

double M, t, deltat;
vector<vector<double>> table;

unsigned seed = (unsigned) std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double max(double a, double b) {
    return (b < a) ? a : b;
}

// returns exponential value of the number in parameter
double exp(double n) {
    return pow(E, n);
}

// u.i.i.d. generator
double get_uniform() {
    std::uniform_real_distribution<double> distribution(0.0, 1.0);
    double number = distribution(generator);
    return (number);
}

/* normal random variate generator */
double box_muller(double m, double s)    // mean m, standard deviation s
{
    double x1, x2, w, y1;
    static double y2;
    static int use_last = 0;

    // use value from previous call
//    if (use_last)
//    {
//        y1 = y2;
//        use_last = 0;
//    }
//    else
    {
        do {
            x1 = 2.0 * get_uniform() - 1.0;
            x2 = 2.0 * get_uniform() - 1.0;
            w = x1 * x1 + x2 * x2;
        } while (w >= 1.0);

        w = sqrt((-2.0 * log(w)) / w);
        y1 = x1 * w;
        y2 = x2 * w;
//        use_last = 1;
    }
    return (m + y1 * s);
}

// normal distribution
double N(const double &z) {
    if (z > 6.0) {
        return 1.0;
    }
    if (z < -6.0) {
        return 0.0;
    }
    double b1 = 0.31938153;
    double b2 = -0.356563782;
    double b3 = 1.781477937;
    double b4 = -1.821255978;
    double b5 = 1.330274429;
    double p = 0.2316419;
    double c2 = 0.3989423;
    double a = fabs(z);
    double t = 1.0 / (1.0 + a * p);
    double b = c2 * exp((-z) * (z / 2.0));
    double n = ((((b5 * t + b4) * t + b3) * t + b2) * t + b1) * t;
    n = 1.0 - b * n;
    if (z < 0.0)
        n = 1.0 - n;
    return n;
}


double normalCDF(double x, double mu, double stddev) {
    return 0.5 * (1 + erf((x - mu) / (stddev * sqrt(2))));
}

double poisson(double theta) {
    double p = exp(-theta);
    double F = p;
    double N = 0.0;
//    std::random_device rd;
//    std::default_random_engine gen(rd());
//    std::uniform_real_distribution<double> u(0.0, 1.0);

    while (get_uniform() > F) {
        N = N + 1.0;
        p = p * theta / N;
        F = F + p;
    }

    return N;
}

double gamma(double a) { //This generates gamma (a, 1)
//    std::random_device rd;
//    std::default_random_engine generator(rd());
//    std::uniform_real_distribution<double> u(0.0, 1.0);
    if (a > 1) {
        double a_bar = a - 1.0;
        double b = (a - (1.0 / (6 * a))) / a_bar;
        double m = 2.0 / a_bar;
        double d = m + 2.0;
        double v = 1.0;
        int accept = -1;

        while (accept < 0) {
            double u1 = get_uniform();
//            double u1 = box_muller(0,1);
//            double u2 = box_muller(0,1);
            double u2 = get_uniform();
            v = b * u2 / u1;
            if (m * u1 - d + v + (1 / v) <= 0.0) { accept = 1; }
            else if (m * log(u1) - log(v) + v - 1.0 <= 0.0) { accept = 1; }
        }
        return a_bar * v;;
    } else {
        double b = (a + exp(1.0)) / exp(1.0);
        int accept = -1;
        double z = 1.0;
        while (accept < 0) {
            double u1 = get_uniform();
//            double u1 = box_muller(0,1);
//            double u2 = box_muller(0,1);
            double u2 = get_uniform();
            double y = b * u1;
            if (y <= 1.0) {
                z = pow(y, 1.0 / a);
                if (u2 <= exp(-z)) { accept = 1; }
            } else {
                z = -log((b - y) / a);
                if (u2 <= pow(z, a - 1.0)) { accept = 1; }
            }
        }

        return z;
    }
}

double chi_square(double v) {
    double a = v / 2.0;
    double beta = 2.0;

    double x = gamma(a);

    return beta * x;
}

//rlist[0] = 0.0;
//vector<double> cir(int numberT,x0, T, alpha, b, sigma)
//{
//
//}

// 0 for the first and rlist[0] - r0 from i = 1 i <numberT


vector<double> exact_euler_CIR(int n, double x0, double T, double alpha, double b, double sigma) { //following pag 124 Glasserman
//    double alpha = k;//This is to preserve the notation used at glasserman using the parameters given at the paper.
    // double b = a / k;
    double d = 4.0 * b * alpha / (sigma * sigma);

//    random_device rd;
//    default_random_engine gen(rd());
//    normal_distribution<double> dist(0.0, 1.0);//N(0,1);

    vector<double> v;
    x0 = sqrt(max(x0, 0));

    v.push_back(x0);
    double delta = T / (double) n;

    if (d > 1.0) {
        for (int j = 0; j < n; j++) {
            double c = sigma * sigma * (1.0 - exp(-alpha * delta)) / (4.0 * alpha);
            double lambda = v.back() * exp(-alpha * delta) / c;
//            double z = get_uniform();
            double z = box_muller(0,1);
            double x = chi_square(d - 1.0);
            v.push_back(c * ((z + sqrt(lambda)) * (z + sqrt(lambda)) + x));
//            cout << v[j] << endl;
        }
    } else {
        for (int j = 0; j < n; j++) {
            double c = sigma * sigma * (1.0 - exp(-alpha * delta)) / (4.0 * alpha);
            double lambda = v.back() * exp(-alpha * delta) / c;
            double N = poisson(lambda / 2.0);
            double x = chi_square(d + 2.0 * N);
            v.push_back(c * x);
//            cout << v[j] << endl;
        }
    }

    ofstream oFile;
    oFile.open("hw3_3.csv", ios::out | ios::ate);
    oFile << "numberT" << "," << "Rate" << endl;

    for (int i = 0; i < n; i++) {
    oFile << i << "," << v[i] << endl;
//    cout << v[i] << endl;
    }
    oFile.close();

    return v;
}

int main() {

    cout << "Cox-Ingersoll-Ross (CIR) Model " << endl;
    cout << "See .csv file for details " << endl;

    exact_euler_CIR(50, 0.3, 1.0, 0.1, 0.4, 2.0);


    cout << endl << endl;
    system("pause");
}
