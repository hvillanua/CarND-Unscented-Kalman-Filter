#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

const float EPSYLON = 1.0e-3;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.0;
  
  // Process noise covariance matrix
  Q_ = MatrixXd(2,2);
  Q_ << std_a_*std_a_, 0,
        0, std_yawdd_*std_yawdd_;

  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  is_initialized_ = false;

  time_us_ = 0.0;

  n_x_ = x_.size();
  n_aug_ = n_x_ + 2;
  lambda_ = 3 - n_aug_;

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Initialize weights
  weights_ = VectorXd(n_sig_);
  float c1 = lambda_ + n_aug_;
  weights_(0) = lambda_ / c1;
  weights_.tail(2*n_aug_).fill(0.5/c1);

  // Sigma points prediction matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  R_radar_ = MatrixXd(3, 3);
  R_radar_.fill(0.0);
  R_radar_(0,0) = std_radr_ * std_radr_;
  R_radar_(1,1) = std_radphi_ * std_radphi_;
  R_radar_(2,2) = std_radrd_ * std_radrd_;

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_.fill(0.0);
  R_lidar_(0,0) = std_laspx_ * std_laspx_;
  R_lidar_(1,1) = std_laspy_ * std_laspy_;

  NIS_radar_ = 0;
  NIS_lidar_ = 0;

  NISvals_radar_.open( "../NIS/NISvals_radar.txt", ios::out );
  NISvals_laser_.open( "../NIS/NISvals_laser.txt", ios::out );

}

UKF::~UKF() {
  NISvals_radar_.close();
  NISvals_laser_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if ((meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) ||
  (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_)) {

    if (!is_initialized_) {
      // first measurement
      cout << "UKF: " << endl;

      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

        /**
        Convert radar from polar to cartesian coordinates and initialize state.
        */
        float rho = meas_package.raw_measurements_[0];
        float phi = meas_package.raw_measurements_[1];
        float rho_dot = meas_package.raw_measurements_[2];

        // transform from polar to cartesian
        float x = rho * cos(phi);
        x = (x < 0) ? min(x, -EPSYLON) : max(x, EPSYLON);
        float y = rho * sin(phi);
        y = (y < 0) ? min(y, -EPSYLON) : max(y, EPSYLON);
        float vx = rho_dot * cos(phi);
        float vy = rho_dot * sin(phi);
        float v = sqrt(vx * vx + vy * vy);

        x_ << x, y, v, 0, 0;
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
        /**
        Initialize state.
        */
        x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      }

      //state covariance matrix P
      P_ = MatrixXd::Identity(n_x_, n_x_);

      // done initializing, no need to predict or update
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
      cout<<"UKF initialized!"<<endl;
      return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    // calculate time difference
    float dt = (meas_package.timestamp_ - time_us_) / 1000000.0; //dt - expressed in seconds
    time_us_ = meas_package.timestamp_;

    // if readings come with the same timestamp dt will be 0, we can skip the predicting step
    if(dt > 0)
    {
      Prediction(dt);
    }

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Radar updates
      UpdateRadar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      // Laser updates
      UpdateLidar(meas_package);
    }
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*****************************************************************************
   *  Sigma points generation
   ****************************************************************************/

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  //create augmented mean state
  x_aug << x_, 0, 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug.bottomRightCorner(2, 2) = Q_;

  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points

  float prod = sqrt(lambda_ + n_aug_);
  MatrixXd x_mat = MatrixXd(n_aug_, n_aug_);
  x_mat << x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug;
  Xsig_aug.col(0) = x_aug;
  Xsig_aug.middleCols(1, n_aug_) = x_mat + prod * A;
  Xsig_aug.middleCols(1+n_aug_, n_aug_) = x_mat - prod * A;

  /*****************************************************************************
   *  Sigma points prediction
   ****************************************************************************/

  //predict sigma points
  double p_x;
  double p_y;
  double v;
  double yaw;
  double yawd;
  double nu_a;
  double nu_yawdd;

  float v_yaw;
  float yawd_t;
  float nu_dt;

  for (int i=0; i<n_sig_; i++)
  {
    p_x = Xsig_aug(0,i);
    p_y = Xsig_aug(1,i);
    v = Xsig_aug(2,i);
    yaw = Xsig_aug(3,i);
    yawd = Xsig_aug(4,i);
    nu_a = Xsig_aug(5,i);
    nu_yawdd = Xsig_aug(6,i);

    v_yaw = v/yawd;
    yawd_t = yawd * delta_t;
    nu_dt = 0.5 * delta_t * delta_t;

    VectorXd H = VectorXd(n_x_);
    VectorXd nu = VectorXd(n_x_);

    nu << nu_dt * cos(yaw) * nu_a,
          nu_dt * sin(yaw) * nu_a,
          delta_t * nu_a,
          nu_dt * nu_yawdd,
          delta_t * nu_yawdd;

    //avoid division by zero
    if (fabs(yawd) > 0.001)
    {
      H << v_yaw * (sin(yaw + yawd_t) - sin(yaw)),
           v_yaw * (-cos(yaw + yawd_t) + cos(yaw)),
           0,
           yawd_t,
           0;
    }
    else
    {
      H << v * cos(yaw) * delta_t,
           v * sin(yaw) * delta_t,
           0,
           yawd_t,
           0;
    }
    //write predicted sigma points into right column
    Xsig_pred_.col(i) << Xsig_aug.col(i).head(n_x_) + H + nu;
  }

  /*****************************************************************************
   *  Predict mean and covariance
   ****************************************************************************/

  //predict state mean
  x_ = Xsig_pred_ * weights_;

  //predict state covariance matrix
  VectorXd diff = VectorXd(n_x_);
  P_.fill(0.0);
  for (int i=0; i<n_sig_; i++)
  {
    diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    diff(3) = fmod(diff(3), M_PI);
    P_ += weights_(i) * diff * diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{

  //set measurement dimension, lidar can measure px and py
  int n_z = 2;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, n_sig_);
  UpdateUKF(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{

  /*****************************************************************************
   *  Predict measurement
   ****************************************************************************/

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);

  float px;
  float py;
  float v;
  float yaw;
  float yawd;

  float rho;
  float phi;
  float rhod;

  for (int i=0; i<n_sig_; i++)
  {
    px = Xsig_pred_(0,i);
    px = (px < 0) ? min(px, -EPSYLON) : max(px, EPSYLON);
    py = Xsig_pred_(1,i);
    py = (py < 0) ? min(py, -EPSYLON) : max(py, EPSYLON);
    v = Xsig_pred_(2,i);
    yaw = Xsig_pred_(3,i);
    yawd = Xsig_pred_(4,i);

    rho = sqrt(px*px + py*py);
    phi = atan2(py, px);
    rhod = (px*cos(yaw)*v + py*sin(yaw)*v) / rho;

    Zsig.col(i) << rho, phi, rhod;
  }

  UpdateUKF(meas_package, Zsig, n_z);
}

void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z)
{

  /*****************************************************************************
   *  Predict measurement
   ****************************************************************************/

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i<n_sig_; i++)
  {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate innovation covariance matrix S
  S.fill(0.0);
  VectorXd diff = VectorXd(n_z);

  MatrixXd R;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    R = R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    R = R_lidar_;
  }

  for (int i=0; i<n_sig_; i++)
  {
    diff = Zsig.col(i) - z_pred;
    // angle normalization
    diff(1) = fmod(diff(1), M_PI);
    S += weights_(i) * diff * diff.transpose();
  }

   S += R;

   /*****************************************************************************
    *  Update
    ****************************************************************************/

   MatrixXd Tc = MatrixXd(n_x_, n_z);
   VectorXd z = meas_package.raw_measurements_;

   //calculate cross correlation matrix
   VectorXd diff1 = VectorXd(n_x_);
   VectorXd diff2 = VectorXd(n_x_);
   Tc.fill(0.0);
   for (int i=0; i<n_sig_; i++)
   {
     diff1 = Xsig_pred_.col(i) - x_;
     diff2 = Zsig.col(i) - z_pred;

     if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
       diff2(1) = fmod(diff2(1), M_PI);
     }

     diff1(3) = fmod(diff1(3), M_PI);

     Tc += weights_(i) * diff1 * diff2.transpose();
     //angle normalization
     Tc(3) = fmod(Tc(3), M_PI);
   }

   //calculate Kalman gain K;
   MatrixXd K = Tc * S.inverse();
   VectorXd z_diff = z-z_pred;

   if (meas_package.sensor_type_ == MeasurementPackage::RADAR){ // Radar
     z_diff(1) = fmod(z_diff(1), M_PI);
   }

   // calculate NIS
   double NIS = z_diff.transpose() * S.inverse() * z_diff;
   if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
   {
     NIS_radar_ = NIS;
     NISvals_radar_ << NIS_radar_ << endl;
   }
   else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
   {
     NIS_lidar_ = NIS;
     NISvals_laser_ << NIS_lidar_ << endl;
   }

   //update state mean and covariance matrix
   x_ += K * z_diff;
   P_ += -K * S * K.transpose();
}
