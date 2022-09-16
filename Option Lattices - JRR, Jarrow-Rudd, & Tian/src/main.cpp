// IE 525 Computational Homework 1
// Tree Methods and Enhancements for Option Pricing
//
// main.cpp
//
// Created by Joseph Loss on 2/17/2019
//
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <chrono>
#include <vector>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

// formatting for tablular output in Question 2
const char separator = ' ';
const int nameWidth = 10;
const int numWidth = 8;

template<typename T> void printElement(T t, const int& width)
{
    cout << left << setw(width) << setfill(separator) << t;
}

// the variables below are for use in QUESTION 1B only.
// inputs to other methods should be adjusted by changing the function arguments within int main{}
int N = 25;
static double **matrix;
double R, up_factor, uptick_prob;
double BlackScholes, Richardson_Extrapolation;
double stock_price = 100;
double strike_price = 105;
double rf_rate = 0.03;
double div_yield = 0.0;
double volatility = 0.20;
double T = 1.0;


double max(double a, double b) {
    return (b < a) ? a : b;
}

double max2(double a, double b, double c) {
    return ((a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c));
}

double NormDist(const double& z) {               // normal distribution function
    if (z > 6.0) { return 1.0; };                // overflow guard
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
    double b = c2 * exp((-z)*(z / 2.0));
    double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
    n = 1.0 - b * n;
    if (z < 0.0) n = 1.0 - n;
    return n;
}

void reinit(double **matrix)               // clearing and reinitializing the price array
{
    for (int i = 0; i <= N; i++)
        for (int j = 0; j <= 2 * N + 1; j++)
            matrix[i][j] = -1;
}

double get_Black_Scholes_Put_Price(const double& stock_price, const double& strike_price,
                                   const double& rf_rate, const double& volatility, const double& T)
{
    double time_sqrt = sqrt(T);
    double d1 = (log(stock_price / strike_price) + rf_rate * T) / (volatility*time_sqrt) + 0.5*volatility*time_sqrt;
    double d2 = d1 - (volatility*time_sqrt);

    return strike_price * exp(-rf_rate * T)*NormDist(-d2) - stock_price * NormDist(-d1);
}

double get_Black_Scholes_Put_Delta(const double& stock_price, const double& strike_price,
                                   const double& rf_rate, const double& volatility, const double& T)
{
    double time_sqrt = sqrt(T);
    double d1 = (log(stock_price / strike_price) + rf_rate * T) /
                (volatility*time_sqrt) + 0.5*volatility*time_sqrt;
    double delta = -NormDist(-d1);
    return delta;
};

double CRR_Binomial_American_Put(double stock_price, double strike_price, double rf_rate, double div_yield, double volatility, double T, int N)
{
    int i, j;
    double uptick_prob;                                 // probability of up movement
    static double St[2500][2500] = { 0.0 };               // stock price at node i,j
    static double put_price[2500][2500] = { 0.0 };        // put price at node i,j
    double a;
    double num = 0.0;
    double up_move = 0.0; double down_move = 0.0;
    double delta_t = 0.0;

    delta_t = T / N;                                      // size of time step
    up_move = exp(volatility*sqrt(delta_t));
    down_move = 1 / up_move;

    a = exp((rf_rate - div_yield)*delta_t);               // growth rate in prob
    uptick_prob = (a - down_move) / (up_move - down_move);

    // compute stock price at each node and initialize call prices
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            St[i][j] = stock_price * (pow(up_move, j)) * (pow(down_move, i - j));
            put_price[i][j] = 0.0;
        }
    }

    // compute payoffs at expiration
    for (int j = N; j >= 0; j--) {
        put_price[N][j] = max(strike_price - St[N][j], 0.0);
    }

    // work backwards
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i; j >= 0; j--) {
            put_price[i][j] = exp(-rf_rate * delta_t)*(uptick_prob*(put_price[i + 1][j + 1]) +
                                                       (1 - uptick_prob)*(put_price[i + 1][j]));
            put_price[i][j] = max(strike_price - St[i][j], put_price[i][j]);
        }
    }
    return put_price[0][0];
}

double JR_Binomial_American_Put(double stock_price, double strike_price, double rf_rate,
                                double div_yield, double volatility, double T, int N)
{
    int i, j;
    double uptick_prob;                                 // probability of up movement
    static double St[2500][2500] = { 0.0 };               // stock price at node i,j
    static double put_price[2500][2500] = { 0.0 };        // put price at node i,j
    double a;
    double num = 0.0;
    double up_move = 0.0; double down_move = 0.0;
    double delta_t = 0.0;


    delta_t = T / N;                                      // size of time step
    up_move = exp(((rf_rate - ((volatility*volatility) / 2))*delta_t) + volatility * sqrt(delta_t));
    down_move = exp(((rf_rate - ((volatility*volatility) / 2))*delta_t) - volatility * sqrt(delta_t));

    uptick_prob = 0.5;

    // compute stock price at each node and initialize call prices
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            St[i][j] = stock_price * (pow(up_move, j)) * (pow(down_move, i - j));
            put_price[i][j] = 0;
        }
    }

    // compute payoffs at expiration
    for (int j = N; j >= 0; j--) {
        put_price[N][j] = max(strike_price - St[N][j], 0.0);
    }

    // work backwards
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i; j >= 0; j--) {
            put_price[i][j] = exp(-rf_rate * delta_t)*(uptick_prob*(put_price[i + 1][j + 1]) +
                                                       (1 - uptick_prob)*(put_price[i + 1][j]));
            put_price[i][j] = max(strike_price - St[i][j], put_price[i][j]);
        }
    }
    return put_price[0][0];
}

double TIAN_Binomial_American_Put(double initial_stock_price, double strike_price, double rf_rate,
                                  double div_yield, double volatility, double T, int N)
{
    int i, j;
    double uptick_prob;                                 // probability of up movement
    static double St[2500][2500] = { 0.0 };               // stock price at node i,j
    static double put_price[2500][2500] = { 0.0 };        // put price at node i,j
    double a; double v; double nu;
    double num = 0.0;
    double up_move = 0.0; double down_move = 0.0;
    double delta_t = 0.0;

    delta_t = T / N;                                      // size of time step
    a = exp((rf_rate - div_yield)*delta_t);               // growth rate in prob
    nu = exp(volatility*volatility*delta_t);

    up_move = 0.5*a*nu*(nu + 1 + sqrt(nu*nu + 2 * nu - 3));
    down_move = 0.5*a*nu*(nu + 1 - sqrt(nu*nu + 2 * nu - 3));
    uptick_prob = (a - down_move) / (up_move - down_move);

    // compute stock price at each node and initialize call prices
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= i; j++) {
            St[i][j] = initial_stock_price * (pow(up_move, j)) * (pow(down_move, i - j));
            put_price[i][j] = 0;
        }
    }

    // compute payoffs at expiration
    for (int j = N; j >= 0; j--) {
        put_price[N][j] = max(strike_price - St[N][j], 0.0);
    }

    // work backwards
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i; j >= 0; j--) {
            put_price[i][j] = exp(-rf_rate * delta_t)*(uptick_prob*(put_price[i + 1][j + 1]) +
                                                       (1 - uptick_prob)*(put_price[i + 1][j]));
            put_price[i][j] = max(strike_price - St[i][j], put_price[i][j]);
        }
    }
    return put_price[0][0];
}

double CRR_BBS_American_Put(int k, int i, double stock_price, int N)
{
    if (matrix[k][i] != -1)
        return matrix[k][i];
    else if (k == N - 1)
        matrix[k][i] = max((strike_price - stock_price),
                           (get_Black_Scholes_Put_Price(stock_price, strike_price, rf_rate, volatility, T / N)));
    else
        matrix[k][i] = max((strike_price - stock_price),
                           (uptick_prob*CRR_BBS_American_Put(k + 1, i + 1, stock_price*up_factor, N) +
                            (1 - uptick_prob)*CRR_BBS_American_Put(k + 1, i - 1, stock_price / up_factor, N)) / R);
    return matrix[k][i];
}

double Generic_Trinomial_American_Put(double initial_stock_price, double strike_price, double rf_rate,
                                      double div_yield, double volatility, double T, int N)
{
    int i, j;
    double downtick_prob, notick_prob, uptick_prob;
    static double St[2500][2500];                           // stock price at node i, j
    static double put_price[2500][2500];                    // put price at node i, j
    double up_move = 0.0;
    double down_move = 0.0;
    double delta_t = T / N;
    double v = rf_rate - 0.5*volatility*volatility;         // drift

    uptick_prob = 0.33333 + (v / volatility)*sqrt(delta_t / 6);
    downtick_prob = 0.33333 - (v / volatility)*sqrt(delta_t / 6);
    notick_prob = 0.33333;

    up_move = exp(volatility*sqrt(3 * delta_t / 2));
    down_move = 1 / up_move;

    // compute stock prices at every node
    for (i = N; i >= 0; i--)
    {
        for (j = -i; j <= i; j++)
        {
            St[i][j] = initial_stock_price * pow(up_move, j);
        }
    }

    // calculate terminal payoffs
    for (j = N; j >= -N; j--)
    {
        put_price[N][j] = max(strike_price - St[N][j], 0);
    }

    // move backwards through trinomial tree
    for (i = N - 1; i >= 0; i--)
    {
        for (j = i; j >= -i; j--)
        {
            put_price[i][j] = max(exp(-rf_rate * delta_t) * (uptick_prob * put_price[i + 1][j + 1] +
                                                             notick_prob * put_price[i + 1][j] +
                                                             downtick_prob * put_price[i + 1][j - 1]), strike_price - St[i][j]);
        }
    }
    return put_price[0][0];
}


double Three_Dimensional_American_Max(double stock1_price, double stock2_price, double stock3_price,
                                      double strike_price, double rf_rate, double stock1_div,
                                      double stock2_div, double stock3_div,
                                      double stock1_sigma, double stock2_sigma, double stock3_sigma,
                                      double rho12, double rho13, double rho23,
                                      double T, int N)
{
    vector <double> St1, St2, St3, tmp;
    vector <vector<double>> P1vals, tmp2;
    vector <vector <vector<double>>> AM_Max_Option;

    double delta_t = T / N;
    double R = exp(rf_rate*delta_t);
    double mu1 = rf_rate - stock1_div - 0.5*stock1_sigma*stock1_sigma;          // stock1 drift
    double mu2 = rf_rate - stock2_div - 0.5*stock2_sigma*stock2_sigma;          // stock2 drift
    double mu3 = rf_rate - stock3_div - 0.5*stock3_sigma*stock3_sigma;          // stock3 drift
    double dx1 = exp(stock1_sigma * sqrt(delta_t));                             // stock1 state step
    double dx2 = exp(stock2_sigma * sqrt(delta_t));                             // stock2 state step
    double dx3 = exp(stock3_sigma * sqrt(delta_t));                             // stock3 state step

    // probability computations
    double puuu = (1 / R)*0.125*(1 + rho12 + rho13 + rho23 + sqrt(delta_t)*(mu1 / stock1_sigma + mu2 / stock2_sigma + mu3 / stock3_sigma));
    double puud = (1 / R)*0.125*(1 + rho12 - rho13 - rho23 + sqrt(delta_t)*(mu1 / stock1_sigma + mu2 / stock2_sigma - mu3 / stock3_sigma));
    double pudu = (1 / R)*0.125*(1 - rho12 + rho13 - rho23 + sqrt(delta_t)*(mu1 / stock1_sigma - mu2 / stock2_sigma + mu3 / stock3_sigma));
    double pudd = (1 / R)*0.125*(1 - rho12 - rho13 + rho23 + sqrt(delta_t)*(mu1 / stock1_sigma - mu2 / stock2_sigma - mu3 / stock3_sigma));
    double pduu = (1 / R)*0.125*(1 - rho12 - rho13 + rho23 - sqrt(delta_t)*(mu1 / stock1_sigma - mu2 / stock2_sigma - mu3 / stock3_sigma));
    double pdud = (1 / R)*0.125*(1 - rho12 + rho13 - rho23 - sqrt(delta_t)*(mu1 / stock1_sigma - mu2 / stock2_sigma + mu3 / stock3_sigma));
    double pddu = (1 / R)*0.125*(1 + rho12 - rho13 - rho23 - sqrt(delta_t)*(mu1 / stock1_sigma + mu2 / stock2_sigma - mu3 / stock3_sigma));
    double pddd = (1 / R)*0.125*(1 + rho12 + rho13 + rho23 - sqrt(delta_t)*(mu1 / stock1_sigma + mu2 / stock2_sigma + mu3 / stock3_sigma));

    // initialize vectors of stock values
    for (int i = 0; i <= 2 * N + 1; i++)
    {
        St1.push_back(0);
        St2.push_back(0);
        St3.push_back(0);
        tmp.push_back(0);       // temp vector to initialize vector-of-vectors
    }

    for (int i = 0; i <= 2 * N + 1; i++)
    {
        tmp2.push_back(tmp);    // temp vector-of-vectors to initialize v-of-v
    }

    St1[1] = stock1_price / pow(dx1, N);
    St2[1] = stock2_price / pow(dx2, N);
    St3[1] = stock3_price / pow(dx3, N);

    for (int i = 2; i <= 2 * N + 1; i++)
    {
        St1[i] = dx1 * St1[i - 1];
        St2[i] = dx2 * St2[i - 1];
        St3[i] = dx3 * St3[i - 1];
    }

    for (int i = 0; i <= 2 * N + 1; i++)
    {
        P1vals.push_back(tmp);
        AM_Max_Option.push_back(tmp2);
    }

    for (int i = 1; i <= 2 * N + 1; i += 2)
    {
        for (int j = 1; j <= 2 * N + 1; j += 2)
        {
            for (int k = 1; k <= 2 * N + 1; k += 2)
                AM_Max_Option[i][j][k] = max2(max2(St1[i], St2[j], St3[k]) - strike_price, 0.0, 0.0);
        }
    }

    for (int t = 1; t <= N; t++)
    {
        for (int i = t + 1; i <= 2 * N + 1 - t; i += 2)
        {
            for (int j = t + 1; j <= 2 * N + 1 - t; j += 2)
            {
                for (int k = t + 1; k <= 2 * N + 1 - t; k += 2)
                {
                    double value_if_held = puuu * AM_Max_Option[i + 1][j + 1][k + 1] + puud * AM_Max_Option[i + 1][j + 1][k - 1] +
                                           pudu * AM_Max_Option[i + 1][j - 1][k + 1] + pudd * AM_Max_Option[i + 1][j - 1][k - 1] + pduu * AM_Max_Option[i - 1][j + 1][k + 1] +
                                           pdud * AM_Max_Option[i - 1][j + 1][k - 1] + pddu * AM_Max_Option[i - 1][j - 1][k + 1] + pddd * AM_Max_Option[i - 1][j - 1][k - 1];
                    AM_Max_Option[i][j][k] = max2(value_if_held, max2(St1[i], St2[j], St3[k]) - strike_price, 0.0);
                }
            }
        }
    }
    return AM_Max_Option[N + 1][N + 1][N + 1];
}


int main(int argc, char* argv[])
{
    cout << "IE 525 Computational Homework 1 " << endl;
    cout << "Tree Methods and Enhancements for Option Pricing " << endl;
    cout << "Created by Joseph Loss (loss2) on 02/17/2019 " << endl << endl;

    cout << "CRR Binomial American Put Approximation: " << endl;
    // inputs: (stock_price, strike_price, rf_rate, div_yield, volatility, T, N)
    cout << CRR_Binomial_American_Put(stock_price, strike_price, rf_rate, div_yield, volatility, T, N) << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "JR Binomial American Put Approximation: " << endl;
    // inputs: (stock_price, strike_price, rf_rate, div_yield, volatility, T, N)
    cout << JR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, 25) << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "TIAN Binomial American Put Approximation: " << endl;
    // inputs: (stock_price, strike_price, rf_rate, div_yield, volatility, T, N)
    cout << TIAN_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, 25) << endl;
    cout << "=========================================================" << endl << endl;

    cout << "TO MOVE ON TO BINOMIAL APPROXIMATIONS WITH RICHARDSON EXTRAPOLATION.. " << endl;
    system("pause");
    cout << endl;

    // ===============================================================================================================
    // ===============================================================================================================
    // Binomial Approximation Improvements

    N = 300;
    matrix = new double *[N + 1];
    for (int i = 0; i <= N; i++)
    {
        matrix[i] = new double[2 * N + 1];
    }
    double b1[300], b2[300];
    N = 300;

    cout << "Converging values from the Binomial Black Scholes Richardson Extrapolation Method:" << endl;
    for (int i = 4; i <= N; i += 2)
    {
        reinit(matrix);

        double delta_t = T / ((float)i);
        up_factor = exp(volatility * sqrt(delta_t));
        R = exp(rf_rate * delta_t);
        uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
        b1[i] = 2 * CRR_BBS_American_Put(0, i, 100, i);

        delta_t = T / ((float)i / 2);
        up_factor = exp(volatility * sqrt(delta_t));
        R = exp(rf_rate * delta_t);
        uptick_prob = (R - (1 / up_factor)) / (up_factor - (1 / up_factor));
        b1[i] -= CRR_BBS_American_Put(0, i, 100, i / 2);

        cout << "BBSR i = " << i << ": " << b1[i];
        cout << endl;

        //        Richardson_Extrapolation = 2*CRR_BBS_American_Put(0, i, 100, (N-1)) - CRR_BBS_American_Put(0, i, 100, (N-1)/2);
    }
    //    cout << endl<< endl<< Richardson_Extrapolation << endl<<endl;

    cout << "---------------------------------------------------------" << endl;

    cout << "Converging values from the BAMR Method:" << endl;
    for (int i = 4; i <= N; i += 2) {
        reinit(matrix);

        double temp1 = CRR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i / 2 - 1);
        double temp2 = CRR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i / 2);
        double temp3 = CRR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i - 1);
        double temp4 = CRR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i);

        b2[i] = 2.0 * (temp4 / 2.0 + temp3 / 2.0) - (temp1 / 2.0 + temp2 / 2.0);

        cout << "BAMR i = " << i << ": " << b2[i];
        cout << endl;
    }

    cout << "=========================================================" << endl;
    cout << endl << endl;
    cout << "TO MOVE ON TO TRINOMIAL LATTICE AND MULTI-DIMENSIONAL OPTIONS.. " << endl;
    system("pause");
    cout << endl << endl;


    // ===============================================================================================================
    // ===============================================================================================================
    // Trinomial Lattice
    cout << "---------------------------------------------------------" << endl;
    cout << "Trinomial American Put Approximation: " << endl;

    // inputs: (stock_price, strike_price, rf_rate, div_yield, volatility, T, N)
    cout << Generic_Trinomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, 25) << endl;
    cout << "---------------------------------------------------------" << endl << endl;


    // ===============================================================================================================
    // ===============================================================================================================
    // Multidimensional Options

    cout << "Three-Dimensional American Max Option Price (steps=15): " << endl;
    std::clock_t startcputime1 = std::clock();

    // inputs: (s1_price, s2_price, s3_price, rf_rate, s1_dividend, s2_dividend, s3_dividend, s1_vol, s2_vol, s3_vol,
    // s1_s2_correlation, s1_s3_correlation, s2_s3_correlation, T, N)
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 15);
    cout << endl;
    double cpu_duration1 = (std::clock() - startcputime1) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration1 << " seconds [CPU Clock] " << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "Three-Dimensional American Max Option Price (steps=30): " << endl;
    std::clock_t startcputime2 = std::clock();
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 30);
    cout << endl;
    double cpu_duration2 = (std::clock() - startcputime2) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration2 << " seconds [CPU Clock] " << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "Three-Dimensional American Max Option Price (steps=60): " << endl;
    std::clock_t startcputime3 = std::clock();
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 60);
    cout << endl;
    double cpu_duration3 = (std::clock() - startcputime3) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration3 << " seconds [CPU Clock] " << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "Three-Dimensional American Max Option Price (steps=120): " << endl;
    std::clock_t startcputime4 = std::clock();
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 120);
    cout << endl;
    double cpu_duration4 = (std::clock() - startcputime4) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration4 << " seconds [CPU Clock] " << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "Three-Dimensional American Max Option Price (steps=240): " << endl;
    std::clock_t startcputime5 = std::clock();
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 240);
    cout << endl;
    double cpu_duration5 = (std::clock() - startcputime5) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration5 << " seconds [CPU Clock] " << endl;
    cout << "---------------------------------------------------------" << endl;

    cout << "Three-Dimensional American Max Option Price (steps=480): " << endl;
    std::clock_t startcputime6 = std::clock();
    cout << Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, 480);
    cout << endl;
    double cpu_duration6 = (std::clock() - startcputime6) / (double)CLOCKS_PER_SEC;
    cout << "Computation finished in " << cpu_duration6 << " seconds [CPU Clock] " << endl;

    cout << endl;
    cout << "=========================================================" << endl;

    cout << "CPU Duration vs Time-Steps" << endl;
    // tabular output of CPU duration vs time-steps
    printElement("Steps: ", nameWidth);
    printElement("  15", numWidth);
    printElement("  30", numWidth);
    printElement("  60", numWidth);
    printElement("  120", numWidth);
    printElement("  240", numWidth);
    printElement("  480", numWidth);
    cout << endl;
    printElement("CPU Time:  ", nameWidth);
    printElement(cpu_duration1, numWidth);
    printElement(cpu_duration2, numWidth);
    printElement(cpu_duration3, numWidth);
    printElement(cpu_duration4, numWidth);
    printElement(cpu_duration5, numWidth);
    printElement(cpu_duration6, numWidth);
    cout << endl;

    cout << "---------------------------------------------------------" << endl;
    cout << "=========================================================" << endl;

    // ===============================================================================================================
    // ===============================================================================================================
    // Output data to csv files

    ofstream oFile;
    oFile.open("question_1a.csv", ios::out | ios::ate);
    oFile << "Lattice Methods" << "," << "," << "Cox-Ross-Rubinstein" << "," << "Jarrow-Rudd" << "," << "Tian" << endl;
    for (int i = 25; i <= 300; i++)
    {
        double temp1, temp2, temp3;
        temp1 = CRR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i);
        temp2 = JR_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i);
        temp3 = TIAN_Binomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i);

        //        *** to display the put price at each step i, uncomment the code below ***

        //        cout << "CRR Binomial American Put Approximation: i=" << i << endl;
        //        cout << temp1 << endl;
        //        cout << "---------------------------------------------------------" << endl;
        //        cout << "JR Binomial American Put Approximation: i=" << i << endl;
        //        cout << temp2 << endl;
        //        cout << "---------------------------------------------------------" << endl;
        //        cout << "TIAN Binomial American Put Approximation: i=" << i << endl;
        //        cout << temp3 << endl;
        //        cout << "=========================================================" << endl;

        double delta_t = T / ((float)i);
        oFile << i << "," << delta_t << "," << temp1 << "," << temp2 << "," << temp3 << endl;
    }
    oFile.close();

    oFile.open("question_1b.csv", ios::out | ios::ate);
    oFile << "step " << "," << "BBSR" << "," << "BAMR" << endl;
    for (int i = 4; i <= 300; i += 2)
    {
        oFile << i << "," << b1[i] << "," << b2[i] << endl;
    }
    oFile.close();

    oFile.open("question_1c.csv", ios::out | ios::ate);
    oFile << "step " << "," << "Tri" << endl;
    for (int i = 25; i <= 300; i++)
    {
        double temp1;
        temp1 = Generic_Trinomial_American_Put(100, 105, 0.03, 0.0, 0.20, 1.0, i);

        //      *** to display the put price at each step i, uncomment the code below ***

        //        cout << "=========================================================" << endl;
        //        cout << "Trinomial American Put Approximation: i=" << i << endl;
        //        cout << temp1 << endl;
        //        cout << "=========================================================" << endl;

        oFile << i << "," << temp1 << endl;
    }
    oFile.close();

    oFile.open("question_2.csv", ios::out | ios::ate);
    oFile << "step,Max,Time" << endl;
    for (int i = 15; i <= 480; i *= 2)
    {
        //        cout << "Three-Dimensional American Max Option Price (steps=" << i << "): " << endl;
        std::clock_t startcputime = std::clock();
        double temp1 = Three_Dimensional_American_Max(100, 100, 100, 100, 0.05, 0.0, 0.0, 0.0, 0.20, 0.20, 0.20, 0.10, 0.10, 0.10, 0.5, i);
        //        cout << temp1 << endl;
        //        cout << endl;
        double cpu_duration = (std::clock() - startcputime) / (double)CLOCKS_PER_SEC;
        //        cout << "Computation finished in " << cpu_duration << " seconds [CPU Clock] " << endl;
        //        cout << "---------------------------------------------------------" << endl;

        oFile << i << "," << temp1 << "," << cpu_duration << endl;
    }
    oFile.close();


    system("pause");
}

