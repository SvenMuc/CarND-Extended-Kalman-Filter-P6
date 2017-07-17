#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // predict the state and the state covariance
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // update the state by using Kalman Filter equations
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  
  long x_size = x_.size();
  MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // update the state by using Extended Kalman Filter equations
  VectorXd y = z - h(x_);
  
  // normalize phi to +/- pi
  y[1] = atan2(sin(y[1]), cos(y[1]));

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();

  long x_size = x_.size();
  MatrixXd I = Eigen::MatrixXd::Identity(x_size, x_size);
  
  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

Eigen::VectorXd KalmanFilter::h(const Eigen::VectorXd &x) {
  // map radar predicted state vector to measurement vector
  float px = x(0);
  float py = x(1);
  float vx = x(2);
  float vy = x(3);
  
  // h(x')
  float rho = sqrt(px*px + py*py);
  float phi = atan2(py, px);
  float rho_dot = (px*vx + py*vy) / rho;
  
  // std::cout << "rho: " << rho << "m  phi: " << phi << " rad (" << phi * 180 / 3.141 << "°) rho_dot: " << rho_dot << std::endl;
  
  VectorXd hx = VectorXd(3);
  hx << rho, phi, rho_dot;
  
  return hx;
}
