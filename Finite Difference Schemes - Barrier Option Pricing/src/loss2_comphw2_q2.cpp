//  Computation Homework #2
//  UIUC - IE525 - Spring 2019
//
//  Created by Joseph Loss on 3/26/2019
//
// Question 2

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include<vector>
#include <chrono>
#include <string>
#include <algorithm>
#include <newmatap.h>
#include <newmat.h>
#include <newmatio.h>

using namespace std;

double T, S0, K, Sb, r, sigma;
double N, M, dx, dt, xmin, xmax, rho;
double d1, d2, d3, d4, d5, d6, d7, d8, a, b;
vector<vector<double>> vector_of_prices;

void write_data(string file_name, vector<vector<double>> &vector_of_prices) {
    ofstream outfile(file_name);
    outfile << ",";
    for (size_t j = 0; j <= M; j++) {
        outfile << j * dt << ",";
    }
    outfile << endl;
    for (size_t i = 0; i <= N; i++) {
        outfile << (i) * dx + xmin << ",";
        for (size_t j = 0; j <= M; j++) {
            {
                outfile << vector_of_prices[i][j] << ",";
            }
        }
        outfile << endl;
    }
}

double Nb(float z) {
    if (z > 6.0) { return 1.0; }; // this guards against overflow
    if (z < -6.0) { return 0.0; };
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
    if (z < 0.0) n = 1.0 - n;
    return n;
}

double P() {
    a = pow(Sb / S0, (-1.0 + 2.0 * r / sigma / sigma));
    b = pow(Sb / S0, (1.0 + 2.0 * r / sigma / sigma));

    d1 = (log(S0 / K) + (r + sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d2 = (log(S0 / K) + (r - sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d3 = (log(S0 / Sb) + (r + sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d4 = (log(S0 / Sb) + (r - sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d5 = (log(S0 / Sb) - (r - sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d6 = (log(S0 / Sb) - (r + sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d7 = (log(S0 * K / Sb / Sb) - (r - sigma * sigma / 2.0) * T) / (sigma * sqrt(T));
    d8 = (log(S0 * K / Sb / Sb) - (r + sigma * sigma / 2.0) * T) / (sigma * sqrt(T));

    double out = K * exp(-r * T) * (Nb(d4) - Nb(d2) - a * (Nb(d7) - Nb(d5))) - S0 * (Nb(d3) - Nb(d1) - b * (Nb(d8) - Nb(d6)));
    return out;
}

double grid_function(int i, int j) {
    if (vector_of_prices[i][j] != -9999) {
        return vector_of_prices[i][j];
    }
    vector_of_prices[i][j] = (-grid_function(i - 2, j + 1) + (2.0 / rho + 2.0) * grid_function(i - 1, j + 1) - grid_function(i, j + 1)) - grid_function(i - 2, j) -
                  (2.0 / rho - 2.0) * grid_function(i - 1, j);

    return vector_of_prices[i][j];
}

int main(int argc, char *argv[]) {
    __int64 temp;
    __int64 temp2;
    __int64 temp1;

    cout << "Computational Homework 2.2 " << endl << endl;
    if (argc != 1) {
        return 0;
    } else {
        {
            cout << "Enter M: " << endl;
            cin >> M;
            cout << "Enter N: " << endl;
            cin >> N;
            cout << "Enter Expiration Time (month): " << endl;
            cin >> T;
            cout << "Enter Sb: " << endl;
            cin >> Sb;
            cout << "Enter S0: " << endl;
            cin >> S0;
            cout << "Enter K: " << endl;
            cin >> K;
            cout << "Enter sigma: " << endl;
            cin >> sigma;
            cout << "Enter r: " << endl;
            cin >> r;
        }
    }

    // Calculations
    T = T / 12.0;
    dt = T / M;
    xmin = 0;
    xmax = 4 * S0;
    dx = ((xmax - xmin) / N);
    rho = dt / dx / dx;
    double R = exp(-r * dt);

    // check and verify inputs:
    if (dx <= 0) {
        cout << "Delta_x must be larger than 0 ";
        return 0;
    }
    if (dt <= 0) {
        cout << "Delta_t must be larger than 0 ";
        return 0;
    }
    if (T <= 0) {
        cout << "Expiration Time must be larger than 0 ";
        return 0;
    }
    if (dt > T) {
        cout << "Expiration Time must be larger than Delta_t ";
        return 0;
    }
    if (xmax < xmin) {
        cout << "Maximum x value must be larger than Minimum x value ";
        return 0;
    }

    // initialize grid
    vector<double> tmp;
    for (size_t j = 0; j <= M; j++) {
        tmp.push_back(-9999);
    }
    for (size_t i = 0; i <= N; i++) {
        vector_of_prices.push_back(tmp);
    }
    for (size_t j = 0; j <= M; j++) {
        vector_of_prices[0][j] = 0.0;
        vector_of_prices[1][j] = 0.0;
        vector_of_prices[N][j] = 0.0;
    }
    for (size_t i = 0; i <= N; i++) {
        if ((i) * dx <= Sb) {
            vector_of_prices[i][M] = 0.0;
        } else if ((i) * dx >= K) {
            vector_of_prices[i][M] = 0;
        } else {
            vector_of_prices[i][M] = K - ((i) * dx);
        }
    }

    Matrix C(N, N), D(N, N);
    ColumnVector g0(N), g1(N), U1(N), U0(N);

    for (size_t i = 1; i <= N; i++) {
        for (size_t j = 1; j <= N; j++) {
            if (i == j) {
                C(i, j) = 1.0 - (r + sigma * sigma * i * i) / 2.0 * dt;
                D(i, j) = 1.0 + (r + sigma * sigma * i * i) / 2.0 * dt;
            } else if ((i - j) == 1) {
                C(i, j) = -(r * i - sigma * sigma * i * i) / 4.0 * dt;
                D(i, j) = (r * i - sigma * sigma * i * i) / 4.0 * dt;
            } else if ((i - j) == -1) {
                C(i, j) = (r * i + sigma * sigma * i * i) / 4.0 * dt;
                D(i, j) = -(r * i + sigma * sigma * i * i) / 4.0 * dt;
            } else {
                C(i, j) = 0;
                D(i, j) = 0;
            }
        }
    }
    for (size_t i = 1; i <= N; i++) {
        g1(i) = 0;
        g0(i) = 0;
        U1(i) = vector_of_prices[i][M];
        U0(i) = 0;
    }
    for (int j = M - 1; j >= 0; j--) {
        g1(1) = vector_of_prices[1][j + 1];
        g1(N) = vector_of_prices[N][j + 1];
        g0(1) = vector_of_prices[1][j];
        g0(N) = vector_of_prices[N][j];
        U0 = D.i() * (C * U1 + rho * (g1 + g0));
        for (size_t i = 1; i <= N; i++) {
            if (i * dx <= Sb) {
                vector_of_prices[i][j] = 0;
            } else {
                if (U0(i) < 0) {
                    vector_of_prices[i][j] = 0;
                } else
                    vector_of_prices[i][j] = U0(i);
            }
            U1(i) = vector_of_prices[i][j];
        }
    }


    cout << endl << endl;
    cout << "DAO Put Price using Crank-Nicolson scheme = " << vector_of_prices[(S0) / dx][0] << endl;

    write_data("vector_of_prices2.csv", vector_of_prices);
    
    cout << "DAO Put Price using exact solution = " << P() << endl;
    cout << "---------------------------------------------" << endl;
    cout << "Difference = " << (vector_of_prices[(S0) / dx][0]) - (P()) << endl;
    cout << endl << endl;
    //cout << P() << endl;

    system("pause");
}