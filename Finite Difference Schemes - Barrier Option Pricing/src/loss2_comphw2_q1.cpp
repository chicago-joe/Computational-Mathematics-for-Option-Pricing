//  Computation Homework #2
//  UIUC - IE525 - Spring 2019
//
//  Created by Joseph Loss on 3/26/2019
//
// Question 1

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <time.h>
#include <vector>
#include <chrono>
#include <string>
#include <algorithm>

using namespace std;

double T, dt, xmin, xmax, dx, c, M, N;
vector<vector<double>> table;

void write_data(string file_name, vector<vector<double>> &table) {
    ofstream outfile(file_name);
    outfile << ",";
    for (size_t j = 0; j <= M; j++) {
        outfile << j * dt << ",";
    }
    outfile << endl;

    for (size_t i = M; i <= N; i++) {
        outfile << (i) * dx + xmin - M * dx << ",";
        for (size_t j = 0; j <= M; j++) {
            if (table[i][j] != -9999) {
                outfile << table[i][j] << ",";
            }
        }
        outfile << endl;
    }
}

double grid_function(int i, int j) {
    if (table[i][j] != -9999) {
        return table[i][j];
    }
    table[i][j] = c * dt / dx * (grid_function(i - 1, j - 1) - grid_function(i, j - 1)) + grid_function(i, j - 1);
    return table[i][j];
}

int main(int argc, char *argv[]) {
//    __int64 temp;
//    __int64 temp1;
//    __int64 temp2;

    cout << "Computational Homework 2.1 " << endl << endl;
    if (argc != 1) {
        return 0;
    } else {
        cout << "Enter Expiration Time (years): " << endl;
        cin >> T;
        cout << "Enter Delta_t: " << endl;
        cin >> dt;
        cout << "Enter Delta_x: " << endl;
        cin >> dx;
        cout << "Enter Minimum x value: " << endl;
        cin >> xmin;
        cout << "Enter Maximum x value: " << endl;
        cin >> xmax;
        cout << "Enter c: " << endl;
        cin >> c;
    }

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
        cout << "Expiration Time (years) must be larger than 0 ";
        return 0;
    }
    if (dt > T) {
        cout << "Expiration Time must be larger than Delta_t ";
        return 0;
    }
    if (dx > 1) {
        cout << "Delta_x must be <= 1 ";
        return 0;
    }
    if (xmax < xmin) {
        cout << "Maximum x value must be larger than Minimum x value ";
        return 0;
    }
    if (ceil((0 - (-1)) / dx) != (0 - (-1)) / dx) {
        cout << "ERROR: Invalid Delta_x value ";
    }
    if (ceil(T / dt) != T / dt) {
        cout << "ERROR: Invalid Delta_t value ";
    }

    // initialize grid
    vector<double> tmp;
    M = ceil(T / dt);
    N = ceil((xmax - (xmin)) / dx) + ceil(T / dt);
    for (size_t j = 0; j <= M; j++) {
        tmp.push_back(-9999);
    }
    for (size_t i = 0; i <= N; i++) {
        table.push_back(tmp);
    }

    for (size_t i = 0; i <= N; i++) {
        if (i * dx + xmin - M * dx < -1) {
            table[i][0] = 0.0;
        } else if (i * dx + xmin - M * dx > 0) {
            table[i][0] = 1.0;
        } else {
            table[i][0] = (double) (i - M) * dx + xmin + 1.0;
        }
    }

    for (size_t i = M; i <= N; i++) {
        grid_function(i, M);
    }

    write_data("question1.csv", table);
    cout << endl << endl;
    system("pause");
}