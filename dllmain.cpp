// dllmain.cpp : Defines the entry point for the DLL application.
#define EXPORTING_DLL 
#include "header.h"
#include "utils.h"
#include <random>

BOOL APIENTRY DllMain( HMODULE hModule,
                       DWORD  ul_reason_for_call,
                       LPVOID lpReserved
                     )
{
    switch (ul_reason_for_call)
    {
    case DLL_PROCESS_ATTACH:
    case DLL_THREAD_ATTACH:
    case DLL_THREAD_DETACH:
    case DLL_PROCESS_DETACH:
        break;
    }
    return TRUE;
}

vector<vector<double>> NormalVariablesMatrix(int N, int Nsim) {

    default_random_engine generator;
    normal_distribution<double> distribution(0.0, 1.0);

    vector<vector<double>> Z;



    for (int i{}; i < Nsim; i++) {

        vector<double> Zi;

        for (int j{ 0 }; j < N; j++) {

            Zi.push_back(distribution(generator));
        }

        Z.push_back(Zi);
    }

    return Z;
}

double CallLookBack(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    double payoff{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS{ S0 }, S{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S, r, vol, h, Zij);

            MinS = min(MinS, S);

        }

        updatePayoff(payoff, S, MinS, Nsim);

    }

    return price(r, T, payoff);
}

double VegaCallLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * vol;
    const double vol1 = vol + eta, vol2 = vol - eta;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS1{ S0 }, MinS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol1, h, Zij);
            updateUnderlying(S2, r, vol2, h, Zij);

            MinS1 = min(MinS1, S1);
            MinS2 = min(MinS2, S2);

        }

        updatePayoff(payoff1, S1, MinS1, Nsim);
        updatePayoff(payoff2, S2, MinS2, Nsim);

    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double vega = priceDerivative(price1, price2, eta);
    return vega;
}


double DeltaCallLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * S0;

    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS1 =  S0 + eta, MinS2 = S0 - eta;
        double S1 =  S0 + eta, S2 = S0 - eta;

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h, Zij);
            updateUnderlying(S2, r, vol, h, Zij);

            MinS1 = min(MinS1, S1);
            MinS2 = min(MinS2, S2);

        }

        updatePayoff(payoff1, S1, MinS1, Nsim);
        updatePayoff(payoff2, S2, MinS2, Nsim);

    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double delta = priceDerivative(price1, price2, eta);
    return delta;
}

double RhoCallLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * r;
    const double r1 = r + eta, r2 = r - eta;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS1{ S0 }, MinS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r1, vol, h, Zij);
            updateUnderlying(S2, r2, vol, h, Zij);

            MinS1 = min(MinS1, S1);
            MinS2 = min(MinS2, S2);

        }

        updatePayoff(payoff1, S1, MinS1, Nsim);
        updatePayoff(payoff2, S2, MinS2, Nsim);

    }

    double price1 = price(r1, T, payoff1);
    double price2 = price(r2, T, payoff2);
    double rho = priceDerivative(price1, price2, eta);
    return rho;
}


double ThetaCallLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double eta = 0.001 * T;
    const double T1 = T + eta, T2 = T - eta;
    const double h1 = T1 / N, h2 = T2 / N;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS1{ S0 }, MinS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h1, Zij);
            updateUnderlying(S2, r, vol, h2, Zij);

            MinS1 = min(MinS1, S1);
            MinS2 = min(MinS2, S2);

        }

        updatePayoff(payoff1, S1, MinS1, Nsim);
        updatePayoff(payoff2, S2, MinS2, Nsim);

    }

    double price1 = price(r, T1, payoff1);
    double price2 = price(r, T2, payoff2);
    double theta = priceDerivative(price1, price2, eta);
    return theta;
}


double GammaCallLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double eta = 0.001 * S0;
    const double h = T / N;
    double payoff1{ 0 }, payoff2{ 0 }, payoff3{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MinS1 = S0 + eta, MinS2 = S0, MinS3 = S0 - eta;
        double S1 = S0 + eta, S2 = S0, S3 = S0 - eta;

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h, Zij);
            updateUnderlying(S2, r, vol, h, Zij);
            updateUnderlying(S3, r, vol, h, Zij);

            MinS1 = min(MinS1, S1);
            MinS2 = min(MinS2, S2);
            MinS3 = min(MinS3, S3);

        }

        updatePayoff(payoff1, S1, MinS1, Nsim);
        updatePayoff(payoff2, S2, MinS2, Nsim);
        updatePayoff(payoff3, S3, MinS3, Nsim);
    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double price3 = price(r, T, payoff3);

    double gamma = priceSecondDerivative(price1, price2, price3, eta);
    return gamma;
}

double PutLookBack(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    double payoff{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS{ S0 }, S{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S, r, vol, h, Zij);

            MaxS = max(MaxS, S);

        }

        updatePayoff(payoff, MaxS, S, Nsim);

    }

    return price(r, T, payoff);
}

double VegaPutLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * vol;
    const double vol1 = vol + eta, vol2 = vol - eta;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS1{ S0 }, MaxS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol1, h, Zij);
            updateUnderlying(S2, r, vol2, h, Zij);

            MaxS1 = max(S1, MaxS1);
            MaxS2 = max(S2, MaxS2);

        }

        updatePayoff(payoff1, MaxS1, S1, Nsim);
        updatePayoff(payoff2, MaxS2, S2, Nsim);

    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double vega = priceDerivative(price1, price2, eta);
    return vega;
}


double DeltaPutLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * S0;

    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS1 = S0 + eta, MaxS2 = S0 - eta;
        double S1 = S0 + eta, S2 = S0 - eta;

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h, Zij);
            updateUnderlying(S2, r, vol, h, Zij);

            MaxS1 = max(S1, MaxS1);
            MaxS2 = max(S2, MaxS2);

        }

        updatePayoff(payoff1, MaxS1, S1, Nsim);
        updatePayoff(payoff2, MaxS2, S2, Nsim);

    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double delta = priceDerivative(price1, price2, eta);
    return delta;
}

double RhoPutLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double h = T / N;
    const double eta = 0.001 * r;
    const double r1 = r + eta, r2 = r - eta;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS1{ S0 }, MaxS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r1, vol, h, Zij);
            updateUnderlying(S2, r2, vol, h, Zij);

            MaxS1 = max(S1, MaxS1);
            MaxS2 = max(S2, MaxS2);

        }

        updatePayoff(payoff1, MaxS1, S1, Nsim);
        updatePayoff(payoff2, MaxS2, S2, Nsim);

    }

    double price1 = price(r1, T, payoff1);
    double price2 = price(r2, T, payoff2);
    double rho = priceDerivative(price1, price2, eta);
    return rho;
}


double ThetaPutLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double eta = 0.001 * T;
    const double T1 = T + eta, T2 = T - eta;
    const double h1 = T1 / N, h2 = T2 / N;
    double payoff1{ 0 }, payoff2{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS1{ S0 }, MaxS2{ S0 };
        double S1{ S0 }, S2{ S0 };

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h1, Zij);
            updateUnderlying(S2, r, vol, h2, Zij);

            MaxS1 = max(S1, MaxS1);
            MaxS2 = max(S2, MaxS2);

        }

        updatePayoff(payoff1, MaxS1, S1, Nsim);
        updatePayoff(payoff2, MaxS2, S2, Nsim);

    }

    double price1 = price(r, T1, payoff1);
    double price2 = price(r, T2, payoff2);
    double theta = priceDerivative(price1, price2, eta);
    return theta;
}


double GammaPutLookBackDF(double S0, double r, double vol, double T, int N, int Nsim, vector<vector<double>> Z)
{
    const double eta = 0.001 * S0;
    const double h = T / N;
    double payoff1{ 0 }, payoff2{ 0 }, payoff3{ 0 };


    for (int i{ 0 }; i < Nsim; i++)
    {
        double MaxS1 = S0 + eta, MaxS2 = S0, MaxS3 = S0 - eta;
        double S1 = S0 + eta, S2 = S0, S3 = S0 - eta;

        for (int j{ 0 }; j < N; j++)
        {

            double Zij = Z[i][j];

            updateUnderlying(S1, r, vol, h, Zij);
            updateUnderlying(S2, r, vol, h, Zij);
            updateUnderlying(S3, r, vol, h, Zij);

            MaxS1 = max(S1, MaxS1);
            MaxS2 = max(S2, MaxS2);
            MaxS3 = max(S3, MaxS3);

        }

        updatePayoff(payoff1, MaxS1, S1, Nsim);
        updatePayoff(payoff2, MaxS2, S2, Nsim);
        updatePayoff(payoff3, MaxS3, S3, Nsim);
    }

    double price1 = price(r, T, payoff1);
    double price2 = price(r, T, payoff2);
    double price3 = price(r, T, payoff3);

    double gamma = priceSecondDerivative(price1, price2, price3, eta);
    return gamma;
}


int SendCallResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim) {

    double callprice;
    double deltacall;
    double vegacall;
    double rhocall; 
    double thetacall; 
    double gammacall;
    int npoints;


    vector<vector<double>> Z = NormalVariablesMatrix(N, Nsim);
        
    callprice = CallLookBack(S0, r,vol, T, N, Nsim, Z);
    deltacall = DeltaCallLookBackDF(S0, r, vol, T, N, Nsim, Z);
    vegacall = VegaCallLookBackDF(S0, r, vol, T, N, Nsim, Z);
    rhocall = RhoCallLookBackDF(S0, r, vol, T, N, Nsim, Z);
    thetacall = ThetaCallLookBackDF(S0, r, vol, T, N, Nsim, Z);
    gammacall = GammaCallLookBackDF(S0, r, vol, T, N, Nsim, Z);

    results[0] = callprice;
    results[1] = deltacall;
    results[2] = vegacall;
    results[3] = rhocall;
    results[4] = thetacall;
    results[5] = gammacall;

    npoints = (int)(3*S0 / 20);

    for (int i{ 0 }; i <= npoints; i++) {
        points[i] = (S0 / 2) + (double) 10 * i;
        pricepoints[i] = CallLookBack(points[i], r, vol, T, N, Nsim, Z);
        deltapoints[i] = DeltaCallLookBackDF(points[i], r, vol, T, N, Nsim, Z);
    }

    return 0;

}



int SendPutResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim) {

    double putprice;
    double deltaput;
    double vegaput;
    double rhoput;
    double thetaput;
    double gammaput;
    int npoints;


    vector<vector<double>> Z = NormalVariablesMatrix(N, Nsim);

    putprice = PutLookBack(S0, r, vol, T, N, Nsim, Z);
    deltaput = DeltaPutLookBackDF(S0, r, vol, T, N, Nsim, Z);
    vegaput = VegaPutLookBackDF(S0, r, vol, T, N, Nsim, Z);
    rhoput = RhoPutLookBackDF(S0, r, vol, T, N, Nsim, Z);
    thetaput = ThetaPutLookBackDF(S0, r, vol, T, N, Nsim, Z);
    gammaput = GammaPutLookBackDF(S0, r, vol, T, N, Nsim, Z);

    results[0] = putprice;
    results[1] = deltaput;
    results[2] = vegaput;
    results[3] = rhoput;
    results[4] = thetaput;
    results[5] = gammaput;

    npoints = (int)(3*S0 / 20);

    for (int i{ 0 }; i <= npoints; i++) {
        points[i] = (S0 / 2) + (double) 10 * i;
        pricepoints[i] = PutLookBack(points[i], r, vol, T, N, Nsim, Z);
        deltapoints[i] = DeltaPutLookBackDF(points[i], r, vol, T, N, Nsim, Z);
    }

    return 0;

}

