#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
using namespace std;

void updateUnderlying(double& S, double r, double vol, double h, double Zij);

void updatePayoff(double& payoff, double S1, double S2, int Nsim);

double price(double r, double T, double payoff);

double priceDerivative(double price1, double price2, double eta);

double priceSecondDerivative(double price1, double price2, double price3, double eta);

