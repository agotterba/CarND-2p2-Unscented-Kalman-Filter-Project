#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  //
  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  //initialize in init
  VectorXd x_;
  VectorXd x_aug_; //[pos1 pos2 vel_abs yaw_angle yaw_rate nu_acceleration nu_yaw_rate]

  ///* state covariance matrix
  //initialize in init
  MatrixXd P_;
  MatrixXd P_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  MatrixXd Xsig_aug_;

  // translated sigma points matrix, prediction, and covariance
  MatrixXd Zsig_;
  VectorXd z_pred_;
  
  //matricies for calculating update; declared separately so they aren't constantly re-allocated
  MatrixXd S_radar_;
  MatrixXd S_laser_;

  MatrixXd R_radar_;
  MatrixXd R_laser_;

  ///* time when the state is true, in us
  //initialize in init
  long long previous_timestamp_;
  long long current_timestamp_;
  //elapsted time from previous measurment, in seconds
  double dt_;
  //count measurements for debugging, and keeping track of percentage of predictions that exceed nis limit
  long long meas_num_;
  long long meas_num_radar_;
  long long meas_num_laser_;
  //count of how many measurements of each type exceeded nis 5% limit
  long long excess_nis_radar_;
  long long excess_nis_laser_;
  //constants of 5% nis limits; defined in init
  double nis_limit_radar_;
  double nis_limit_laser_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;
  int n_radar_;
  int n_laser_;
  
  ///* Augmented state dimension
  int n_aug_;
  int n_aug_points_;

  ///* Sigma point spreading parameter
  double lambda_;


  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_pack);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction();

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_pack);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_pack);

  //Initalize the kalman filter
  void Initalize(MeasurementPackage meas_pack);

  //normalize angle to between -PI and +PI
  double Normalize_angle (double theta);

  //Generate sigma points
  void Gen_sig_points();
  //Predict values (mean) and covariance for sigma points
  void Predict_sig_points();
  //predict new mean and covariance
  void Predict_xP();
  //Initialize weight values for all sigma points
  void Set_weights();
  //Convert predictions into radar measurement space
  void Predict_radar();
  //Convert predictions into lidar measurement space
  void Predict_lidar();
  //Calculate new x_ values and covariance for new radar measurement
  void Calc_radar(VectorXd z);
  //Calculate new x_ values and covariance for new lidar measurement
  void Calc_lidar(VectorXd z);
  //Find NIS for a radar measurement, and check against limit
  double Nis_radar(VectorXd z);
  //Fine NIS for a lidar measurement, and check against limit
  double Nis_lidar(VectorXd z);
};

#endif /* UKF_H */
