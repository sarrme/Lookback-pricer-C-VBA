#define WIN32_LEAN_AND_MEAN             
#include <windows.h>
#include <string>

#ifdef EXPORTING_DLL

extern "C" __declspec(dllexport) int SendCallResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim);
extern "C" __declspec(dllexport) int SendPutResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim);

#else

extern "C" __declspec(dllimport) int SendCallResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim);
extern "C" __declspec(dllimport) int SendPutResultsToVBA(double* results, double* points, double* pricepoints, double* deltapoints, double S0, double r, double vol, double T, int N, int Nsim);
#endif // EXPORTING_DLL

