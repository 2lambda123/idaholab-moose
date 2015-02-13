/*
 *      This class contains a set of multi-dimensional interpolation functions for both scattered data and data lying on a cartesian grid
 *
 *      Sources:
 *      - General
 *        * Numerical Recipes in C++ 3rd edition
 *      - MD spline
 *        * Christian Habermann, Fabian Kindermann, "Multidimensional Spline Interpolation: Theory and Applications", Computational Economics, Vol.30-2, pp 153-169 (2007) [http://link.springer.com/article/10.1007%2Fs10614-007-9092-4]
 *      - Inverse distance weighting
 *        * http://en.wikipedia.org/wiki/Inverse_distance_weighting
 *
 */


#ifndef ND_INTERPOLATION_FUNCTIONS_H
#define ND_INTERPOLATION_FUNCTIONS_H

#include <vector>
#include <string>
#include <iostream>
#include <stdio.h>

bool checkUpperBound(double upper_bound, std::vector<double> values);
bool checkLowerBound(double lower_bound, std::vector<double> values);
std::vector<double> int2binary(int value, int size);

class NDInterpolation
{
public:
  virtual double interpolateAt(std::vector<double> point_coordinate);
  virtual double getGradientAt(std::vector<double> point_coordinate);
  virtual void   fit(std::vector< std::vector<double> > coordinates, std::vector<double> values);
  std::vector<double> NDinverseFunction(double F_min, double F_max);
  std::vector<double> NDinverseFunctionGrid(double F);

  double averageCellValue(std::vector<double> center, std::vector<double> dx);

  void updateRNGparameters(double tolerance, double initial_divisions);

  double NDderivative(std::vector<double> coordinate);
  double integral(std::vector<double> coordinate, int samples=1000);

  int returnDimensionality(){return _dimensions;};

  NDInterpolation();
  ~NDInterpolation();

protected:
  std::string _data_filename;
  bool _completed_init;
  int _dimensions;

  double _tolerance;
  int _initial_divisions;

  std::vector<double> _cellPoint0;
  std::vector<double> _cellDxs;

  double minkowskiDistance(std::vector<double> point1, std::vector<double> point2, double p);
  double vectorNorm(std::vector<double> point, double p);

  bool pivotCellCheck(std::vector<std::vector<double> > cell, double F);
  int vertexOutcome(std::vector<double> vertex, double F);
  void cellsFilter(std::vector<std::vector<std::vector<double> > >& vertices, double F);
  void refinedCellDivision(std::vector<std::vector<std::vector<double> > >& refinedCell, std::vector<std::vector<double> > cell, int divisions);
  std::vector<int> arrayConverter(int oneDcoordinate, int divisions, int n_dimensions);
  std::vector<std::vector<double> > generateNewCell(std::vector<int> NDcoordinate, std::vector<double> coordinateOfPointZero, std::vector<double> dxs, int n_dimensions);
  std::vector<std::vector<double> > pickNewCell(std::vector<std::vector<std::vector<double> > > cellsSet, int seed);
  std::vector<double> getCellCenter(std::vector<std::vector<double> > cell);

  double OneDderivative(double fxph, double fx, double fxmh);
  double derivativeStep(std::vector<double> coordinate, int loop);
};

class NDSpline: public NDInterpolation
{
public:
  double interpolateAt(std::vector<double> point_coordinate);
  double getGradientAt(std::vector<double> point_coordinate);
  double integralSpline(std::vector<double> point_coordinate);
  double spline_cartesian_marginal_integration(double coordinate,int marginal_variable);
  double spline_cartesian_inverse_marginal(double CDF,int marginal_variable, double precision);

  void   fit(std::vector< std::vector<double> > coordinates, std::vector<double> values);

  NDSpline(std::string filename);
  NDSpline(std::string filename, std::vector<double> alpha, std::vector<double> beta);
  NDSpline(std::vector< std::vector<double> > & discretizations, std::vector<double> & values, std::vector<double> alpha, std::vector<double> beta);

  //std::vector< std::vector<double> > getDiscretizations(){
  //	  std::cout<<"but why!"<< std::endl;
  //	  return _discretizations;};

  void getDiscretizations(std::vector< std::vector<double> > & vector){
	  for(int i=0; i<_discretizations.size();i++){
		  std::vector<double> temp;
		  for(int j=0; j<_discretizations.at(i).size(); j++)
			  temp.push_back(_discretizations.at(i).at(j));
		  vector.push_back(temp);
	  }
	  std::cout<< "xxxx "<< vector.size() << std::endl;
  }

  NDSpline();
  ~NDSpline();

  bool checkUB(double upper_bound);
  bool checkLB(double lower_bound);

  bool checkBoundaries(std::vector<double> point);

private:
  std::vector< std::vector<double> > _discretizations;
  std::vector<double> _values;
  std::vector<double> _spline_coefficients;
  std::vector<double> _hj;
  std::vector<double> _alpha;
  std::vector<double> _beta;

  std::vector<double> _min_disc;
  std::vector<double> _max_disc;

  //void initializeCoefficientsVector();
  void saveCoefficient(double value, std::vector<int> coefficient_coordinate);
  double retrieveCoefficient(std::vector<int> coefficient_coordinate);

  double spline_cartesian_interpolation(std::vector<double> point_coordinate);
  double spline_cartesian_integration(std::vector<double> point_coordinate);
  double getPointAtCoordinate(std::vector<int> coordinates);

  int fromNDto1Dconverter(std::vector<int> coordinate);
  std::vector<int> from1DtoNDconverter(int oneDcoordinate, std::vector<int> indexes);

  void calculateCoefficients();
  std::vector<double> fillArrayCoefficient(int n_dimensions, std::vector<double> & data, std::vector<int> & loop_locator);

  void from2Dto1Drestructuring(std::vector<std::vector<double> > & twoDdata, std::vector<double> & oneDdata);
  void from1Dto2Drestructuring(std::vector<std::vector<double> > & twoDdata, std::vector<double> & oneDdata, int spacing);

  double phi(double t);
  double PHI(double t);
  double u_k(double x, std::vector<double> & discretizations, double k);
  double U_K(double x, std::vector<double> & discretizations, double k);
  void tridag(std::vector<double> & a, std::vector<double> & b, std::vector<double> & c, std::vector<double> & r, std::vector<double> & u);
  std::vector<double> getCoefficients(std::vector<double> & y, double h, double alpha, double beta);
  //void iterationStep(int nDim, std::vector<double> & coefficients, std::vector<double> & data);

  std::vector<double> coefficientRestructuring(std::vector<std::vector<double> > matrix);
  std::vector<std::vector<double> > tensorProductInterpolation(std::vector<std::vector<double> > step1, double h, double alpha, double beta);
  std::vector<std::vector<double> > matrixRestructuring(std::vector<std::vector<double> > step1);
  std::vector<double> getValues(std::vector<int> & loop_locator);

  double val1(double t);
  double val2(double t);
  double val3(double t);
  double val4(double t);
  double val5(double t);
  double val6(double t);

};

class InverseDistanceWeighting: public NDInterpolation
{
public:
  double interpolateAt(std::vector<double> point_coordinate);
  double getGradientAt(std::vector<double> point_coordinate);
  void   fit(std::vector< std::vector<double> > coordinates, std::vector<double> values);
  std::vector<double> NDinverseFunction(double F_min, double F_max);
  InverseDistanceWeighting(std::string filename, double p);
  InverseDistanceWeighting(double p);
  bool checkUB(double upper_bound);
  bool checkLB(double lower_bound);

  std::vector<double> get_cellPoint0(){return _cellPoint0;};
  std::vector<double> get_cellDxs(){return _cellDxs;};

private:
  int _number_of_points;
  double _p;
  std::vector<double> _values;
  //std::vector<double> _cellPoint0;
  //std::vector<double> _cellDxs;
  std::vector< std::vector<double> > _point_coordinates;
};

//class NDlinear: public NDInterpolation
//{
//public:
//  double interpolateAt(std::vector<double> point_coordinate);
//  NDlinear(std::string filename);
//  NDlinear();
//  ~NDlinear();
//  double linear_interpolation(std::vector<double> point_coordinate);
//  std::vector<double> getValues(std::vector<int> & loopLocator);
//  int fromNDto1Dconverter(std::vector<int> coordinate);
//  bool checkBoundaries(std::vector<double> point);
//  std::vector<int> from1DtoNDconverter(int oneDcoordinate, std::vector<int> indexes);
//
//  bool checkUB(double upper_bound);
//  bool checkLB(double lower_bound);
//
//private:
//  int _dimensions;
//  int _number_of_points;
//  bool _completedInit;
//  double _p;
//  std::vector<double> _values;
//  std::vector<double> _minDisc;
//  std::vector<double> _maxDisc;
//  std::vector< std::vector<double> > _discretizations;
//};

class MicroSphere: public NDInterpolation
{
public:
  double interpolateAt(std::vector<double> point_coordinate);
  double getGradientAt(std::vector<double> point_coordinate);
  void   fit(std::vector< std::vector<double> > coordinates, std::vector<double> values);
  MicroSphere(std::string filename, double p, int precision);
  MicroSphere(double p, int precision);
private:
  int _number_of_points;
  double _p;
  std::vector<double> _values;
  std::vector< std::vector<double> > _point_coordinates;
  int _precision;
  std::vector< std::vector<double> > _unit_vector;
  void MSinitialization();
  double cosValueBetweenVectors(std::vector<double> point1, std::vector<double> point2);
};


#endif
