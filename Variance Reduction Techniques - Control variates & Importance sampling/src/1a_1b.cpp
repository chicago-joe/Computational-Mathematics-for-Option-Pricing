// IE525 Computational Homework 4
// UIUC - IE525 - Spring 2019
//
// Created by Joseph Loss on 4/21/2019
//
// Question 1a & 1b

#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
using namespace std;

double r, K, s0, T, sigma, B;
int no_of_simulations, no_of_barriers;

unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator(seed);

double max(double a, double b) {
	return (b < a) ? a : b;
}

// u.i.i.d. generator
double get_uniform()
{
	std::uniform_real_distribution <double> distribution(0.0, 1.0);
	double number = distribution(generator);
	return (number);
}

double N(const double& z) {
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
	double b = c2 * exp((-z)*(z / 2.0));
	double n = ((((b5*t + b4)*t + b3)*t + b2)*t + b1)*t;
	n = 1.0 - b * n;
	if (z < 0.0) n = 1.0 - n;
	return n;
}

double option_price_call_black_scholes(const double& S,       // spot (underlying) price
	const double& K,       // strike (exercise) price,
	const double& r,       // interest rate
	const double& sigma,   // volatility 
	const double& time) {  // time to maturity 
	double time_sqrt = sqrt(time);
	double d1 = (log(S / K) + r * time) / (sigma*time_sqrt) + 0.5*sigma*time_sqrt;
	double d2 = d1 - (sigma*time_sqrt);
	return S * N(d1) - K * exp(-r * time)*N(d2);
}


int main()
{
	cout << "European Down-and-Out Discrete Barrier Options Pricing via Monte Carlo Simulation" << endl << endl;
	
	// PARAMETERS:
	T = 25.0 / 252.0;			
	r = 0.03;
	sigma = 0.60;
	s0 = 99.0;
	K = 105.0;					// strike price
	B = 90.0;					// barrier price
	no_of_simulations = 10000.0;
	no_of_barriers = 25.0;

	cout << "--------------------------------" << endl;
	cout << "PARAMETERS: " << endl;
	cout << "Expiration Time (Years) = " << T << endl;
	cout << "Risk Free Interest Rate = " << r << endl;
	cout << "Volatility (%age of stock value) = " << sigma * 100 << endl;
	cout << "Initial Stock Price = " << s0 << endl;
	cout << "Strike Price = " << K << endl;
	cout << "Barrier Price = " << B << endl;
	cout << "Number of Trials = " << no_of_simulations << endl;
	cout << "Number of Discrete Barriers = " << no_of_barriers << endl;
	cout << "--------------------------------" << endl << endl;

	double *simulated_barrier_call_option_price = new double[no_of_simulations + 1];
	double *squared_simulated_barrier_call_option_price = new double[no_of_simulations + 1];
	double *simulated_vanilla_call_option_price = new double[no_of_simulations + 1];
	double *squared_simulated_vanilla_call_option_price = new double[no_of_simulations + 1];
	double *predicted_barrier_call_option_price = new double[no_of_simulations + 1];
	double cov = 0.0, price_simulated_barrier_call_option = 0.0, variance_simulated_barrier_call_option = 0.0, price_simulated_vanilla_call_option = 0.0, variance_simulated_vanilla_call_option = 0.0;
	double delta_T = T / ((double)no_of_barriers);
	double delta_R = (r - 0.5*pow(sigma, 2))*delta_T;
	double delta_SD = sigma * sqrt(delta_T);

	for (int k = 0; k < no_of_simulations; k++)
	{
		//create 4 paths
		double current_stock_price1 = s0;
		double current_stock_price2 = s0;
		double current_stock_price3 = s0;
		double current_stock_price4 = s0;

		//these will store the barrier-breach-adjusted price
		double S1 = 0, S2 = 0, S3 = 0, S4 = 0;

		//tracking if stock hits barrier 
		bool current_stock_price1_hit_barrier = false;
		bool current_stock_price2_hit_barrier = false;
		bool current_stock_price3_hit_barrier = false;
		bool current_stock_price4_hit_barrier = false;


		// create unit-normal variates using Box-Muller
		for (int j = 0; j < no_of_barriers; j++) {
			double x = get_uniform();
			double y = get_uniform();
			double a = sqrt(-2.0*log(x)) * cos(6.283185307999998*y);
			double b = sqrt(-2.0*log(x)) * sin(6.283185307999998*y);

			// if Stock Price hits the barrier price --> Hit_Barrier = TRUE 
			current_stock_price1 = current_stock_price1 * exp(delta_R + delta_SD * a);
			current_stock_price1_hit_barrier = (current_stock_price1 <= B ? true : current_stock_price1_hit_barrier);

			current_stock_price2 = current_stock_price2 * exp(delta_R - delta_SD * a);
			current_stock_price2_hit_barrier = (current_stock_price2 <= B ? true : current_stock_price2_hit_barrier);

			current_stock_price3 = current_stock_price3 * exp(delta_R + delta_SD * b);
			current_stock_price3_hit_barrier = (current_stock_price3 <= B ? true : current_stock_price3_hit_barrier);

			current_stock_price4 = current_stock_price4 * exp(delta_R - delta_SD * b);
			current_stock_price4_hit_barrier = (current_stock_price4 <= B ? true : current_stock_price4_hit_barrier);
		}

		simulated_vanilla_call_option_price[k] = (max(0.0, current_stock_price1 - K) +
			max(0.0, current_stock_price2 - K) +
			max(0.0, current_stock_price3 - K) +
			max(0.0, current_stock_price4 - K)) / 4.0;

		price_simulated_vanilla_call_option += simulated_vanilla_call_option_price[k]; 

		squared_simulated_vanilla_call_option_price[k] = (pow(max(0.0, current_stock_price1 - K), 2) +
			pow(max(0.0, current_stock_price2 - K), 2) +
			pow(max(0.0, current_stock_price3 - K), 2) +
			pow(max(0.0, current_stock_price4 - K), 2)) / 4.0;

		variance_simulated_vanilla_call_option += squared_simulated_vanilla_call_option_price[k]; 

		// If barrier_hit = TRUE --> Invalidate Price Path (= 0)
		S1 = current_stock_price1_hit_barrier == true ? 0 : current_stock_price1;
		S2 = current_stock_price2_hit_barrier == true ? 0 : current_stock_price2;
		S3 = current_stock_price3_hit_barrier == true ? 0 : current_stock_price3;
		S4 = current_stock_price4_hit_barrier == true ? 0 : current_stock_price4;

		// Compute avg payoff of the simulated price paths that didn't hit the barrier
		simulated_barrier_call_option_price[k] = (max(0.0, S1 - K) +
			max(0.0, S2 - K) +
			max(0.0, S3 - K) +
			max(0.0, S4 - K)) / 4.0;

		price_simulated_barrier_call_option += simulated_barrier_call_option_price[k]; 

		squared_simulated_barrier_call_option_price[k] = (pow(max(0.0, S1 - K), 2) +
			pow(max(0.0, S2 - K), 2) +
			pow(max(0.0, S3 - K), 2) +
			pow(max(0.0, S4 - K), 2)) / 4.0;

		variance_simulated_barrier_call_option += squared_simulated_barrier_call_option_price[k]; 

		cov += (max(0.0, current_stock_price1 - K)*max(0.0, S1 - K) +
			max(0.0, current_stock_price2 - K)*max(0.0, S2 - K) +
			max(0.0, current_stock_price3 - K)*max(0.0, S3 - K) +
			max(0.0, current_stock_price4 - K)*max(0.0, S4 - K)) / 4.0;
	}

	price_simulated_barrier_call_option = exp(-r * T)*(price_simulated_barrier_call_option / ((double)no_of_simulations));		 // discounted payoff of barrier option
	price_simulated_vanilla_call_option = exp(-r * T)*(price_simulated_vanilla_call_option / ((double)no_of_simulations));		 // discounted payoff of vanilla option

	double barrier_SE = sqrt((variance_simulated_barrier_call_option - pow(price_simulated_barrier_call_option, 2)) / (no_of_simulations - 1) / no_of_simulations);
	double vanilla_SE = sqrt((variance_simulated_vanilla_call_option - pow(price_simulated_vanilla_call_option, 2)) / (no_of_simulations - 1) / no_of_simulations);

	double beta = cov / variance_simulated_vanilla_call_option;
	double alpha = price_simulated_barrier_call_option - beta * price_simulated_barrier_call_option;
	double sum_of_squared_residuals = 0.0;

	for (int k = 0; k < no_of_simulations; k++)
	{
		predicted_barrier_call_option_price[k] = alpha + beta * simulated_vanilla_call_option_price[k];
		sum_of_squared_residuals += pow(predicted_barrier_call_option_price[k] - simulated_barrier_call_option_price[k], 2);
	}


	cout << "The average Barrier Call Price via explicit simulation of price paths = " << price_simulated_barrier_call_option << endl;
	cout << "Standard error in Barrier Call Option Price = " << barrier_SE << endl;
	cout << "Actual error in Barrier Call Option Price = " << price_simulated_barrier_call_option << " - 4.647650 = " << price_simulated_barrier_call_option - 4.647650 << endl;
	cout << "--------------------------------" << endl << endl;
	cout << "The average Vanilla Call Price via explicit simulation of price paths = " << price_simulated_vanilla_call_option << endl;
	cout << "Standard error in Vanilla Call Option Price = " << vanilla_SE << endl;
	cout << "Actual error in Vanilla Call Option Price = " << price_simulated_vanilla_call_option << " - " << option_price_call_black_scholes(s0, K, r, sigma, T) << 
		" = " << price_simulated_vanilla_call_option - option_price_call_black_scholes(s0, K, r, sigma, T) << endl << endl;

	cout << "R-squared = " << 1 - sum_of_squared_residuals / variance_simulated_barrier_call_option << endl;
	cout << "--------------------------------" << endl;

	cout << endl;
	system("pause");

}
