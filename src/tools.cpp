#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0)
  {
    cout<<"Estimations vector has size 0"<<endl;
    return rmse;
  }
  else if(estimations.size() != ground_truth.size())
  {
    cout<<"Estimations vector has different size than ground truth"<<endl;
    return rmse;
  }

  //accumulate squared residuals
  int acc_res = 0;
  VectorXd aux_vec;
  for(int i=0; i < estimations.size(); ++i){
    aux_vec = estimations[i] - ground_truth[i];
    aux_vec = aux_vec.array() * aux_vec.array();
    rmse += aux_vec;
  }

  //calculate the mean
  rmse = rmse/estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}
