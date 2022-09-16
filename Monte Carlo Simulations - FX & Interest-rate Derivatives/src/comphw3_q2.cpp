//  Computational Homework #3
//  UIUC - IE525 - Spring 2019
//
//  Created by Joseph Loss on 4/10/2019
//
// Question 2

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

double spread_monte_carlo_EU_option(double stock1, double stock2, double rf_rate, double sigma1, double sigma2,
                                    double rho, double T, char type, int no_of_trials, int no_of_steps) {
    int i, j;
    double mu1 = (rf_rate - 0.5 * sigma1 * sigma1);
    double mu2 = (rf_rate - 0.5 * sigma2 * sigma2);
    double srho = sqrt(1 - rho * rho);
    double sum1 = 0.0;
    double sum2 = 0.0;
    double St1 = 0.0;
    double St2 = 0.0;
    double deviate1 = 0.0;
    double deviate2 = 0.0;
    double z1 = 0.0;
    double z2 = 0.0;
    double CT = 0.0;
    double SD = 0.0;
    double SE = 0.0;
    double option_value = 0.0;
    double deltat = 0.0;

    //no_of_steps = 1;

    cout << "--------------------------------------------" << endl;
    cout << "Monte-Carlo Exchange Option Pricing Model " << endl;
    cout << "Input Parameters: " << endl;
    cout << "Stock1 initial price (v0) = " << stock1 << endl;
    cout << "Stock2 initial price (u0) = " << stock2 << endl;
    cout << "Risk-free rate = " << rf_rate << endl;
    cout << "Stock1 sigma = " << sigma1 << endl;
    cout << "Stock2 sigma = " << sigma2 << endl;
    cout << "Rho = " << rho << endl;
    cout << "Time to expiration = " << T << endl;
    cout << "Call (C) or Put (P) = " << type << endl;
    cout << "Number of simulations = " << no_of_trials << endl;
    cout << "Number of steps = " << no_of_steps << endl << endl;


    for (int i = 0; i < no_of_trials; i++) {
        // initialize prices for each simulation
        St1 = stock1;
        St2 = stock2;

        for (j = 0; j < no_of_steps; j++) {
            deltat = T / no_of_steps;

            // generate deviates
            deviate1 = box_muller(0.0, 1.0);
            deviate2 = box_muller(0.0, 1.0);

            // calc correlated deviates
            z1 = deviate1;
            z2 = rho * deviate1 + srho * deviate2;
            St1 = St1 * exp(mu1 * deltat + sigma1 * z1 * sqrt(deltat));
            St2 = St2 * exp(mu2 * deltat + sigma2 * z2 * sqrt(deltat));
        }

        if (type == 'C')
            CT = max(St1 - St2, 0);
        // else
        // CT = max(St1 + St2, 0);

        sum1 = sum1 + CT;
        sum2 = sum2 + CT * CT;
    }

    option_value = exp(-rf_rate * T) * (sum1 / no_of_trials);

    SD = sqrt((sum2 - sum1 * sum1 / no_of_trials) * exp(-2 * rf_rate * T) / (no_of_trials - 1));
    SE = SD / sqrt(no_of_trials);

    cout << "Monte-Carlo Exchange Option Price = $ " << option_value << endl;

    ofstream oFile;
    oFile.open("hw3_2.csv", ios::out | ios::ate);
    oFile << "Number of Steps" << "," << "Option Value" << endl;
    for (int l = 1; l <= 1000; l++) {
        double option_val = 0;
        sum1 = 0.0;
            for (int i = 0; i <= l; i++) {
        // initialize prices for each simulation
        St1 = stock1;
        St2 = stock2;

        for (j = 0; j < no_of_steps; j++) {
            deltat = T / no_of_steps;

            // generate deviates
            deviate1 = box_muller(0.0, 1.0);
            deviate2 = box_muller(0.0, 1.0);

            // calc correlated deviates
            z1 = deviate1;
            z2 = rho * deviate1 + srho * deviate2;
            St1 = St1 * exp(mu1 * deltat + sigma1 * z1 * sqrt(deltat));
            St2 = St2 * exp(mu2 * deltat + sigma2 * z2 * sqrt(deltat));
        }

        if (type == 'C')
            CT = max(St1 - St2, 0);
        // else
        // CT = max(St1 + St2, 0);

        sum1 = sum1 + CT;
    }
    option_value = exp(-rf_rate * T) * (sum1 / (l+1));

    oFile << l << "," << option_value << endl;
//    cout << v[i] << endl;
    }
    oFile.close();


    return option_value;
}

double spread_option_black_scholes(double stock1, double stock2, double rf_rate,
        double sigma1, double sigma2, double rho, double T, char type) {

    double sum1 = 0.0;
//    double St1, St2 = 0.0;
    double CT = 0.0;
    double option_value = 0.0;

    double sigmahat = sqrt(sigma1 * sigma1 + sigma2 * sigma2 - 2 * rho * sigma1 * sigma2);
    double d1 = (log(stock1/ stock2) + sigmahat * sigmahat * (T / 2)) / (sigmahat * sqrt(T));
    double d2 = d1 - sigmahat * sqrt(T);

    CT = stock1 * N(d1) + stock2 * -N(d2);
    sum1 = sum1 + CT;

    option_value = exp(-rf_rate * T) * (sum1);
    cout << "Black-Scholes Exchange Option Price = $ " << option_value << endl;



    return option_value;
}

//void write_data(string file_name, vector<vector<double>> &table)
//{
//    ofstream outfile(file_name);
//    outfile << ",";
//    //j to T
//    for (size_t j = 0; j <= no_of_trials; j++)
//    {
//        outfile << j*deltat <<",";
//    }
//    outfile << endl;
//    for (size_t i = 0; i <= N; i++)
//    {
//        //i to S
//        outfile <<(i)*deltax+xmin<< ",";
//        for (size_t j = 0; j <= M; j++)
//        {
//            //if (table[i][j] != -9999)
//            {
//                outfile << table[i][j] << ",";
//            }
//        }
//        outfile << endl;
//    }
//}


int main() {


    spread_monte_carlo_EU_option(50, 60, 0.05, 0.3, 0.4, 0.7, 5.0 / 12.0, 'C', 100000, 100);
    spread_option_black_scholes(50, 60, 0.05, 0.3, 0.4, 0.7, 5.0 / 12.0, 'C');



    cout << endl;
    system("pause");

}
