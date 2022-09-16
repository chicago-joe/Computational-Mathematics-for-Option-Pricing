//  Computation Homework #2
//  UIUC - IE525 - Spring 2019
//
//  Created by Joseph Loss on 3/26/2019
//
// Question 3

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <chrono>
#include <string>
#include <algorithm>
#include <cmath>

using namespace std;

double T, xmin, xmax, M, N, eps, omega, domega, K, r, sigma, k;
vector<vector<double>> vector_of_prices, vector_of_payoffs;
double oldloops, rho, dt, dx;

void write_data(string file_name, vector<vector<double>> &vector_of_prices) {
    ofstream outfile(file_name);
    outfile << ",";
    for (size_t j = 0; j <= M; j++) {
        outfile << -(j * dt / sigma / sigma * 2.0 - T) << ",";
    }
    outfile << endl;
    for (size_t i = 0; i <= N; i++) {
        outfile << K * exp((i) * dx + xmin) << ",";
        for (size_t j = 0; j <= M; j++) {
            {
                outfile << vector_of_prices[i][j] << ",";
            }
        }
        outfile << endl;
    }
}

double pay_off(double x, double tau) {
    return exp(0.5 * (k + 1.0) * (k + 1.0) * (tau)) * max(exp(0.5 * (k - 1.0) * x) - exp(0.5 * (k + 1.0) * x), 0.0);
}

int main(int argc, char *argv[]) {
    cout << "Computational Homework 2.3 " << endl << endl;
    if (argc != 1) {
        return 0;
    } else {
        {
            cout << "Enter Expiration Time (months): " << endl;
            cin >> T;
            cout << "Enter step M: " << endl;
            cin >> M;
            cout << "Enter step N: " << endl;
            cin >> N;
            cout << "Enter eps: " << endl;
            cin >> eps;
            cout << "Enter omega: " << endl;
            cin >> omega;
            cout << "Enter domega: " << endl;
            cin >> domega;
            cout << "Enter Strike Price (K): " << endl;
            cin >> K;
            cout << "Enter Risk-free rate (r): " << endl;
            cin >> r;
            cout << "Enter Sigma: " << endl;
            cin >> sigma;
        }
    }

    T = T / 12;
    oldloops = 10000;
    xmax = log(K * 3.0 / K);
    xmin = log(K / 3.0 / K);
    dx = (xmax - xmin) / N;
    dt = T * sigma * sigma / 2.0 / M;
    rho = dt / dx / dx;
    k = 2.0 * r / sigma / sigma;

    // initialize grids:
    vector<double> x, u, g, b;
    double loop, error, y;

    vector<double> tmp;
    for (size_t j = 0; j <= M; j++) {
        tmp.push_back(-9999);
    }
    for (size_t i = 0; i <= N; i++) {
        vector_of_prices.push_back(tmp);
    }
    for (size_t i = 0; i <= N; i++) {
        tmp.clear();
        for (size_t j = 0; j <= M; j++) {
            tmp.push_back(pay_off(i * dx + xmin, j * dt));
        }
        vector_of_payoffs.push_back(tmp);
    }
    for (size_t j = 0; j <= M; j++) {
        vector_of_prices[0][j] = K * exp(-0.5 * (2.0 * r / sigma / sigma - 1.0) * (0 * dx + xmin) -
                              0.25 * (2.0 * r / sigma / sigma + 1.0) * (2.0 * r / sigma / sigma + 1.0) * j * dt) * vector_of_payoffs[0][j];

        vector_of_prices[N][j] = K * exp(-0.5 * (2.0 * r / sigma / sigma - 1.0) * (N * dx + xmin) -
                              0.25 * (2.0 * r / sigma / sigma + 1.0) * (2.0 * r / sigma / sigma + 1.0) * j * dt) * vector_of_payoffs[N][j];
    }
    for (size_t i = 0; i <= N; i++) {
            vector_of_prices[i][0] = K * exp(-0.5 * (2.0 * r / sigma / sigma - 1.0) * (i * dx + xmin) -
                                  0.25 * (2.0 * r / sigma / sigma + 1.0) * (2.0 * r / sigma / sigma + 1.0) * 0.0 * dt) * vector_of_payoffs[i][0];
    }
    for (size_t i = 0; i <= N; i++) {
        x.push_back(xmin + i * dx);
        u.push_back(vector_of_payoffs[i][0]);
        g.push_back(vector_of_payoffs[i][0]);
        b.push_back(-9999);
    }

    for (int j = 1; j <= M; j++) {
        for (size_t i = 1; i <= N - 1; i++) {
            g[i] = vector_of_payoffs[i][j];
            b[i] = u[i] + rho / 2.0 * (u[i + 1] - 2.0 * u[i] + u[i - 1]);
        }
        g[0] = vector_of_payoffs[0][j];
        g[N] = vector_of_payoffs[N][j];
        u[0] = g[0];
        u[N] = g[N];
        loop = 0;
        do {
            error = 0;
            for (size_t i = 1; i < N; i++) {
                y = (b[i] + rho / 2.0 * (u[i - 1] + u[i + 1])) / (1.0 + rho);
                if (true) {
                    y = max(g[i], u[i] + omega * (y - u[i]));
                } else {
                    y = (u[i] - y) * (u[i] - y);
                }
                error += (u[i] - y) * (u[i] - y);
                u[i] = y;
            }
            loop++;
            if (loop > 10000) {
                cout << "error";
            }
        } while (error > eps);
        if (loop > oldloops) {
            domega *= -1.0;
        }
        
        omega += domega;
        oldloops = loop;
        for (size_t i = 1; i <= N; i++) {
            vector_of_prices[i][j] = K * exp(-0.5 * (k - 1) * (i * dx + xmin) - 0.25 * (k + 1) * (k + 1) * j * dt) * u[i];
            }
    }

    cout << endl << "---------------------------------- " << endl;
    cout << "American Put Price using PSOR method = " << vector_of_prices[ceil((log(K / K) - xmin) / dx)][M] << endl << endl;
    write_data("question3.csv", vector_of_prices);

    system("pause");
}

