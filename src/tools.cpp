#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  Eigen::VectorXd rmse = VectorXd(4);
  rmse << 0,0,0,0;
  if (estimations.size() == 0){
    cout<<"NOTE: CalculateRMSE: called with estimations list of size 0\n";
    return rmse;
  }
  if (ground_truth.size() != estimations.size()){
    cout<<"ERROR: CalculateRMSE: estimations vector has size "<<estimations.size()
        <<" but ground_truth vector has size " <<ground_truth.size()<<" \n";
    return rmse;
  }

  //VectorXd diff = estimations.array() - ground_truth.array();
  //VectorXd diff2 = diff.array() * diff.array();
  //accumulate squared residuals

  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
  //calculate the mean
  // ... your code here
  VectorXd mean;
  rmse = rmse / estimations.size();
  
  //calculate the squared root
  // ... your code here
  rmse = rmse.array().sqrt();
  
  return rmse;
}
