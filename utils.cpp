// dllmain.cpp : Defines the entry point for the DLL application.
#include "utils.h"

void updateUnderlying(double& S, double r, double vol, double h, double Zij) {
    S = S * exp((r - 0.5 * pow(vol, 2)) * h + vol * sqrt(h) * Zij);
}

void updatePayoff(double& payoff, double S1, double S2, int Nsim) {
    payoff += (S1 - S2) / Nsim;
}

double price(double r, double T, double payoff) {
    return exp(-r * T) * payoff;
}

double priceDerivative(double price1, double price2, double eta) {
    return (price1 - price2) / (2 * eta);
}

double priceSecondDerivative(double price1, double price2, double price3, double eta)
{
    return (price1 + price3 - 2 * price2) / pow(eta, 2);
}


