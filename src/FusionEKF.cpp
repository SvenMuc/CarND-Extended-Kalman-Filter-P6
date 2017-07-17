#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;
  previous_timestamp_ = 0;
  
  // measurement matrix - lase
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  // measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;
  
  // measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;
  
  // set acceleration noise
  noise_ax = 9;
  noise_ay = 9;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    // state transition matrix
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ << 1, 0, 1, 0,
               0, 1, 0, 1,
               0, 0, 1, 0,
               0, 0, 0, 1;

    // state covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
               0, 1, 0, 0,
               0, 0, 1000, 0,
               0, 0, 0, 1000;
    /**
     * Measurement format:
     * LASER: [meas_px, meas_py, timestamp, gt_px, gt_py, gt_vx, gt_vy]
     * RADAR: [meas_rho, meas_phi, meas_rho_dot, timestamp, gt_px, gt_py, gt_vx, gt_vy]
     */
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // initialize radar states x = [px, py, vx, vy]
      // convert radar from polar to cartesian coordinates and initialize state.
      float px = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
      float py = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
     
      ekf_.x_<< px, py, 0.0, 0.0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // initialize laser states x = [px, py, vx, vy]
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0.0, 0.0;
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }
  
  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  
  // update state transition matrix F with actual dt
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;
  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  // update process noise covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << (pow(dt, 4.f) / 4.f * noise_ax), 0, (pow(dt, 3.f) / 2.f * noise_ax), 0,
              0, (pow(dt, 4.f) / 4.f * noise_ay), 0, (pow(dt, 3.f) / 2.f * noise_ay),
              (pow(dt, 3.f) / 2.f * noise_ax), 0, (pow(dt, 2.f) * noise_ax), 0,
              0, (pow(dt, 3.f) / 2.f * noise_ay), 0, (pow(dt, 2.f) * noise_ay);
  
  ekf_.Predict();
  
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    
    if(!Hj_.isZero()) {
      ekf_.H_ = Hj_;
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
  } else {
    // laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
  
  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
