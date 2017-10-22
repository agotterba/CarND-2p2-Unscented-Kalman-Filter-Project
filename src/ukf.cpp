#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //vector length  for native, radar, laser, and augmented points
  n_x_ = 5;
  n_radar_ = 3;
  n_laser_ = 2;
  n_aug_ = 7;
  //number of sigma points
  n_aug_points_ = 2 * n_aug_ + 1;
  //lambda
  lambda_ = 3 - n_aug_;
  
  // initial state vector [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  x_ = VectorXd(n_x_);
  x_aug_ = VectorXd(n_aug_);

  // initial covariance matrix 
  P_ = MatrixXd(n_x_, n_x_);
  P_aug_ = MatrixXd(n_aug_, n_aug_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  //std_a_ = 30;
  std_a_ = 2;
  //DEBUG: match coursework
  //std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  //std_yawdd_ = 30;
  std_yawdd_ = 0.75;
  //DEBUG: match coursework
  //std_yawdd_ = 0.2;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Laser nis limit for 5% excess on a measurement with 2 degrees of freedom
  nis_limit_laser_ = 5.991;
  
  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;
  //DEBUG: match coursework
  //std_radphi_ = 0.0175;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DEBUG: match coursework
  //std_radrd_ = 0.1;

  // Radar nis limit for 5% excess on a measurement with 3 degrees of freedom
  nis_limit_radar_ = 7.815;
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  //detects first measurement
  is_initialized_ = false;
  //keep track of measurement numbers for debug and nis percentage excess calculation
  meas_num_ = 0;
  meas_num_radar_ = 0;
  meas_num_laser_ = 0;
  excess_nis_radar_ = 0;
  excess_nis_laser_ = 0;

  //declare these as part of ukf object (and separate for laser and radar),
  //so they aren't being constantly re-allocated
  S_radar_ = MatrixXd(n_radar_,n_radar_);
  S_laser_ = MatrixXd(n_laser_,n_laser_);
  R_radar_ = MatrixXd(n_radar_,n_radar_);
  R_laser_ = MatrixXd(n_laser_,n_laser_);
  
  //define these as members of the object so memory isn't constantly being reallocated
  Xsig_pred_ = MatrixXd(n_x_, n_aug_points_);
  Xsig_aug_  = MatrixXd(n_aug_, n_aug_points_);

  //weights for the different sigma points
  weights_ = VectorXd(n_aug_points_);
  //only set weights once, since they never change, depend only on constants defined above
  Set_weights();
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  current_timestamp_ = meas_pack.timestamp_;

  //DEBUG
  //current_timestamp_ = previous_timestamp_ + 100000; //set to 0.1s to match course work

  dt_ = (current_timestamp_ - previous_timestamp_)/1e6;
  meas_num_ += 1;
  
  if (! is_initialized_){
    Initalize(meas_pack);
    previous_timestamp_ = current_timestamp_;
    //cout <<"Initialized with meas_num_ "<<meas_num_<<"\nx_ = \n"<<x_<<"\nP_ = \n"<<P_<<"\n";
    return;
  }
  //cout<<"Processing meas_num_ "<<meas_num_<<" with dt_ = "<<dt_<<"\n";
  Prediction();
  //cout<<"  After prediction, state is:\nx_=\n"<<x_<<"\nP_=\n"<<P_<<"\n";
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_ ){
      //cout<<"  Updating RADAR measurement\n";
      UpdateRadar(meas_pack);
    }else{
      cout <<"  Skipping RADAR measurement as use_radar_ is false\n";
    }
  }else{ //LIDAR
    if(use_laser_ ){
      //cout<<"  Updating LIDAR measurement\n";
      UpdateLidar(meas_pack);
    }else{
      cout <<"  Skipping LIDAR measurement as use_laser_ is false\n";
    }
  }
  cout<<"  After measurement "<<meas_num_<<", state is:\nx_=\n"<<x_<<"\nP_=\n"<<P_<<"\n";
  previous_timestamp_ = current_timestamp_;
  return;
}

void UKF::Initalize(MeasurementPackage meas_pack){
  VectorXd xynew = VectorXd(2);
  if (meas_pack.sensor_type_ == MeasurementPackage::RADAR){
    double ro = meas_pack.raw_measurements_(0);
    double theta = Normalize_angle(meas_pack.raw_measurements_(1));
    //double ro_dot = meas_pack.raw_measurements_(2);
    xynew(0) = ro * cos(theta);
    xynew(1) = ro * sin(theta);

  }else{
    xynew(0) = meas_pack.raw_measurements_(0);
    xynew(1) = meas_pack.raw_measurements_(1);

  }

   x_ << xynew[0],xynew[1],0,0,0;
   P_ << 0.1,   0,   0,   0,   0,
           0, 0.1,   0,   0,   0,
           0,   0,   1,   0,   0,
           0,   0,   0,   1,   0,
           0,   0,   0,   0,   1;
   //DEBUG to match with lab work
   // x_ <<  5.7441,
   //        1.3800,
   //        2.2049,
   //        0.5015,
   //        0.3528;
   // P_ <<  0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
   //       -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
   //        0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
   //       -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
   //       -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
  is_initialized_ = true;
  return;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction() {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  Gen_sig_points();
  Predict_sig_points();
  Predict_xP();
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  meas_num_laser_ += 1;
  //cout<<"  converting prediction to lidar space\n";
  Predict_lidar();
  double nis = Nis_lidar(meas_pack.raw_measurements_);
  double nis_percent = 100 * (meas_num_laser_ - excess_nis_laser_)/meas_num_laser_;
  //cout<<"nis laser is "<<nis<<"\n";
  cout<<"laser nis percent under limit: "<<nis_percent<<"\n";
  //cout<<"  updating based on lidar measurement\n";
  Calc_lidar(meas_pack.raw_measurements_);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_pack) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  meas_num_radar_ += 1;
  //cout<<"  converting prediction to radar space\n";
  Predict_radar();
  double nis = Nis_radar(meas_pack.raw_measurements_);
  double nis_percent = 100 * (meas_num_radar_ - excess_nis_radar_)/meas_num_radar_;
  //cout<<"nis radar is "<<nis<<"\n";
  cout<<"radar nis percent under limit: "<<nis_percent<<"\n";
  //cout<<"  updating based on radar measurement\n";
  Calc_radar(meas_pack.raw_measurements_);

  // VectorXd debug = VectorXd(3);
  // debug<<
  //   5.9214,   //rho in m
  //   0.2187,   //phi in rad
  //   2.0062;   //rho_dot in m/s
  //Calc_radar(debug);

}

double UKF::Normalize_angle (double theta){
  int count = 0;

  theta = atan2(sin(theta),cos(theta));
  return theta;
  
  // while (theta < -1 * M_PI){
  //   cout << "NOTE: normalize_angle:  adding 2pi to theta of "<<theta<<"\n";
  //   theta += 2* M_PI;
  //   cout <<"  theta is now "<<theta<<"\n";
  //   count = count + 1;
  //   //if (count > 4){
  //   //  cout<<"ERROR:Normalize_angle: iterated too many times\n";
  //   //  return theta;
  //   //}
  // }
  // while (theta > M_PI){
  //   cout << "NOTE: normalize_angle: subtracting 2pi from theta of "<<theta<<"\n";
  //   theta -= 2 * M_PI;
  //   cout <<"  theta is now "<<theta<<"\n";
  //   count = count + 1;
  //   //if (count > 4){
  //   //  cout<<"ERROR:Normalize_angle: iterated too many times\n";
  //   //  return theta;
  //   //}
  // }
}

void UKF::Gen_sig_points(){
  ////cout<<"  generating sigma points\n";
  ////cout<<"    update augmented mean state\n";
  x_aug_.head(n_x_) = x_;
  x_aug_(5) = 0; //acceleration noise has mean of 0
  x_aug_(6) = 0; //yaw noise has mean of 0
  ////cout<<"    built x_aug:\n"<<x_aug_<<"\n";

  ////cout<<"    update augmented covariance matrix\n";
  P_aug_.fill(0.0);
  P_aug_.topLeftCorner(n_x_,n_x_) = P_;
  P_aug_(5,5) = std_a_ * std_a_;
  P_aug_(6,6) = std_yawdd_ * std_yawdd_;
  ////cout<<"    built P_aug:\n"<<P_aug_<<"\n";

  ////cout<<"    update square root matrix\n";
  MatrixXd Plam_aug = P_aug_ * (lambda_ + n_aug_);
  ////cout<<"    built Plam_aug:\n"<<Plam_aug<<"\n";
  MatrixXd sqrtPlam_aug = Plam_aug.llt().matrixL();
  ////cout<<"    built sqrtPl;am_aug:\n"<<sqrtPlam_aug<<"\n";
  
  ////cout<<"    update augmented sigma points\n";
  Xsig_aug_.col(0) = x_aug_;
  for (int i = 0; i < n_aug_; i++){
    Xsig_aug_.col(i+1) = x_aug_ + sqrtPlam_aug.col(i);
    Xsig_aug_.col(n_aug_ + i+1)   = x_aug_ - sqrtPlam_aug.col(i);
  }
  //cout<<"    generated sigma points:\n"<<Xsig_aug_<<"\n";
  return;
}

void UKF::Predict_sig_points(){
  ////cout<<"  Predicting sigma points\n";
  
  double dt2 = 0.5 * dt_ * dt_;
  ////cout<<"    predicting updated value for each sigma point\n";
    
  for (int i = 0; i < n_aug_points_; i++){
    double px      = Xsig_aug_(0,i);
    double py      = Xsig_aug_(1,i);
    double v       = Xsig_aug_(2,i);
    double phi     = Xsig_aug_(3,i);
    double phi_dot = Xsig_aug_(4,i);
    double nu_a    = Xsig_aug_(5,i);
    double nu_phi  = Xsig_aug_(6,i);
    if ((phi_dot > 1e-6)||(phi_dot < -1e-6)){ //phi_dot is not too close to zero
      Xsig_pred_(0,i) = px  + (v / phi_dot) * (sin(phi + phi_dot * dt_) - sin(phi)) + dt2 * cos(phi) * nu_a;   
      Xsig_pred_(1,i) = py  + (v / phi_dot) * (cos(phi) - cos(phi + phi_dot * dt_)) + dt2 * sin(phi) * nu_a;
    }else{
      //cout<<"      phi_dot is close to 0; using alternate equations\n";
      Xsig_pred_(0,i) = px  + v * cos(phi) * dt_ + dt2 * cos(phi) * nu_a;  
      Xsig_pred_(1,i) = py  + v * sin(phi) * dt_ + dt2 * sin(phi) * nu_a;
    }
    double phi_mod = 0;
    Xsig_pred_(2,i) = v       + 0           + dt_ * nu_a;
    //DO NOT NORMALIZE: NEED TO BE > PI FOR INITIAL STEPS OF PREDICTION (DIFF) TO WORK
    Xsig_pred_(3,i) = phi     + phi_dot*dt_ + dt2 * nu_phi + phi_mod;
    Xsig_pred_(4,i) = phi_dot + 0           + dt_ * nu_phi;
  }  
  //cout<<"    Predicted sigma points: Xsig_pred:\n"<<Xsig_pred_<<"\n";
  return;
}

void UKF::Predict_xP(){
  ////cout<<"  Predicting x and P\n";
  x_.fill(0.0);
  P_.fill(0.0);

  ////cout<<"    predict state mean\n";
  for (int i = 0; i < n_aug_points_; i++){
    x_ = x_ + (weights_(i) * Xsig_pred_.col(i));
  }

  ////cout<<"    predict state covariance matrix\n";
  ////cout<<"    x is \n"<<x_<<"\n";
  for (int i = 0; i < n_aug_points_; i++){
    //cout<<"    sig_pred is \n"<<Xsig_pred_.col(i)<<"\n";
    VectorXd diff = Xsig_pred_.col(i) - x_;
    //cout<<"    raw diff is \n"<<diff<<"\n";
    diff(3) = Normalize_angle(diff(3));
    //cout<<"    norm diff is \n"<<diff<<"\n";
    P_ = P_ + (weights_(i) * diff * diff.transpose());
  }
  ////cout<<"    predicted xP.  x is:\n"<<x_<<"\n    and P is:\n"<<P_<<"\n";
  return;
}

void UKF::Set_weights(){
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < n_aug_points_; i++){
    weights_(i) = 0.5/(lambda_ + n_aug_);
  }
  return;
}

void UKF::Predict_radar(){
  Zsig_ = MatrixXd(n_radar_, n_aug_points_);
  ////cout<<"    transforming sigma points into measurement space\n";
  for (int i = 0; i < n_aug_points_; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double phi = Xsig_pred_(3,i);
    double phi_dot = Xsig_pred_(4,i);
    double v1 = cos(phi)*v;
    double v2 = sin(phi)*v;
    Zsig_(0,i) = sqrt(px*px + py*py);
    Zsig_(1,i) = atan2(py,px);
    if (Zsig_(0,i) > 1e-6){
      Zsig_(2,i) = (px*v1 + py*v2) / Zsig_(0,i);
    }else{
      Zsig_(2,i) = 0;
    }
  }
  ////cout<<"    calculating mean predicted measurement\n";
  z_pred_ = VectorXd(n_radar_);
  z_pred_.fill(0.0);
  for (int i = 0; i < n_aug_points_; i++){
      z_pred_ = z_pred_ + weights_(i)*Zsig_.col(i);
  }
  ////cout<<"    calculating covariance matrix S\n";
  S_radar_.fill(0.0);
  for (int i = 0; i < n_aug_points_; i++){
      VectorXd diff = Zsig_.col(i) - z_pred_;
      diff(1) = Normalize_angle(diff(1));
      S_radar_ = S_radar_ + weights_(i) * diff * diff.transpose();
  }
  R_radar_.fill(0.0);
  R_radar_(0,0) = std_radr_   * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_  * std_radrd_;
  S_radar_ = S_radar_ + R_radar_;

  ////cout<<"    Finished converting predictions to radar.  z_pred_ is:\n"<<z_pred_
  ////    <<"\n  Zsig_ is :\n"<<Zsig_<<"\n"
  ////    <<"\n  S_ is:\n"<<S_<<"\n";
  return;
}

void UKF::Calc_radar(VectorXd z){
  ////cout <<"    calculating cross correlation matrix\n";
  ////cout <<"    z is \n"<<z<<"\n";
  MatrixXd Tc = MatrixXd(n_x_,n_radar_);
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_points_; i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = Normalize_angle(x_diff(3));
      VectorXd z_diff = Zsig_.col(i) - z_pred_;
      z_diff(1) = Normalize_angle(z_diff(1));
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  ////cout<<"    calculating Kalman gain K\n";
  MatrixXd Kgain = Tc * S_radar_.inverse();
  ////cout<<"    updating state mean and covariance matrix\n";
  z(1) = Normalize_angle(z(1));
  VectorXd f_diff = z - z_pred_;
  f_diff(1) = Normalize_angle(f_diff(1));
  x_ = x_ + Kgain * f_diff;
  x_(3) = Normalize_angle(x_(3));
  P_ = P_ - Kgain * S_radar_ * Kgain.transpose();

  ////cout<<"    calculated radar update.  x is:\n"<<x_<<"\n    and P is\n"<<P_<<"\n";
}


void UKF::Predict_lidar(){
  Zsig_ = MatrixXd(n_laser_, n_aug_points_);
  ////cout<<"    transforming sigma points into measurement space\n";
  for (int i = 0; i < n_aug_points_; i++){
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    Zsig_(0,i) = px;
    Zsig_(1,i) = py;
  }
  ////cout<<"    calculating mean predicted measurement\n";
  z_pred_ = VectorXd(n_laser_);
  z_pred_.fill(0.0);
  for (int i = 0; i < n_aug_points_; i++){
      z_pred_ = z_pred_ + weights_(i)*Zsig_.col(i);
  }
  ////cout<<"    calculating covariance matrix S\n";
  S_laser_.fill(0.0);
  R_laser_.fill(0.0);
  R_laser_(0,0) = std_laspx_   * std_laspx_;
  R_laser_(1,1) = std_laspy_ * std_laspy_;
  for (int i = 0; i < n_aug_points_; i++){
      VectorXd diff = Zsig_.col(i) - z_pred_;
      S_laser_ = S_laser_ + weights_(i) * diff * diff.transpose();
  }
  S_laser_ = S_laser_ + R_laser_;

  ////cout<<"    Finished converting predictions to laser.  z_pred_ is:\n"<<z_pred_
  ////    <<"\n  Zsig_ is :\n"<<Zsig_<<"\n"
  ////    <<"\n  S_ is:\n"<<S_<<"\n";
  return;
}

void UKF::Calc_lidar(VectorXd z){
  ////cout <<"    calculating cross correlation matrix\n";
  ////cout <<"    z is \n"<<z<<"\n";
  MatrixXd Tc = MatrixXd(n_x_,n_laser_);
  Tc.fill(0.0);
  for (int i = 0; i < n_aug_points_; i++){
      VectorXd x_diff = Xsig_pred_.col(i) - x_;
      x_diff(3) = Normalize_angle(x_diff(3));
      VectorXd z_diff = Zsig_.col(i) - z_pred_;
      Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  ////cout<<"    calculating Kalman gain K\n";
  MatrixXd Kgain = Tc * S_laser_.inverse();
  ////cout<<"    updating state mean and covariance matrix\n";
  x_ = x_ + Kgain * (z - z_pred_);
  x_(3) = Normalize_angle(x_(3));
  P_ = P_ - Kgain * S_laser_ * Kgain.transpose();

  ////cout<<"    calculated laser update.  x is:\n"<<x_<<"\n    and P is\n"<<P_<<"\n";
}

double UKF::Nis_radar(VectorXd z){
  VectorXd z_diff = z - z_pred_;
  double epsilon = z_diff.transpose() * S_radar_.inverse() * z_diff;
  if (epsilon > nis_limit_radar_){
    excess_nis_radar_ += 1;
  }
  return epsilon;
}

double UKF::Nis_lidar(VectorXd z){
  VectorXd z_diff = z - z_pred_;
  double epsilon = z_diff.transpose() * S_laser_.inverse() * z_diff;
  if (epsilon > nis_limit_laser_){
    excess_nis_laser_ += 1;
  }
  return epsilon;
}
