#ifndef DISTRIBUTIONFUNCTIONS_H
#define DISTRIBUTIONFUNCTIONS_H

#include <sstream>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <string>
#include <iostream>
#include <stdio.h>
#include <cmath>   // to use erfc error function
#include <ctime>   // for rand() and srand()

/*
 *  distributionFunctions
 *      source: Numerical Recipes in C++ 3rd edition
 */

//typedef boost::numeric::ublas::matrix<double> matrixDouble;

void matrixConversion(std::vector<std::vector<double> > original, double converted[]);
void matrixBackConversion(double original[], std::vector<std::vector<double> > converted);
void inverseMatrix(double* A, int N);
void computeInverse(const std::vector<std::vector<double> > & matrix, std::vector<std::vector<double> > & inverse);
// convert a matrix stored in a vector to a matrix stored in a vector of vector
void vectorToMatrix(int &rows, int & columns, std::vector<double> &vecMatrix, std::vector<std::vector<double> > &_cov_matrix);
double getDeterminant(std::vector<std::vector<double> > matrix);

void nrerror(const char error_text[]);

double gammp(double a, double x);

double loggam(double xx);


#define ITMAX 100
#define EPSW 3.0e-7

void gser(double *gamser,double a,double x,double *gln);

#define FPMIN 1.0e-30

void gcf(double *gammcf,double a,double x,double *gln);

double gammaFunc(double x);

double betaFunc(double alpha, double beta);

double logGamma(double x);

double betaContFrac(double a, double b, double x);

double betaInc(double a, double b, double x);

double normRNG(double mu, double sigma, double RNG);

void LoadData(double** data, int dimensionality, int cardinality, std::string filename);

double calculateCustomPdf(double position, double fitting, double** dataSet, int numberSamples);

double calculateCustomCDF(double position, double fitting, double** dataSet, int numberSamples);

double rk_gauss();

double STDgammaRNG(double shape);

double gammaRNG(double shape, double scale);

double betaRNG(double alpha, double beta);

double ModifiedLogFunction(double x);

double AbramStegunApproximation(double t);


#endif /* DISTRIBUTIONFUNCTIONS_H */
