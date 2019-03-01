#include "ukf.h"
#include "Eigen/Dense"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2.;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 1.;
  
  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

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

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  is_initialized_ = false;

  /**
   * Prediction step init
   */
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  Q_ = (MatrixXd(2, 2) << std_a_*std_a_, 0, 0, std_yawdd_*std_yawdd_).finished();

  // create vector for weights
  double w = 1/(lambda_+n_aug_);
  weights_ = (VectorXd(2 * n_aug_ + 1) << w*lambda_, w/2.0, w/2.0, w/2.0, w/2.0, w/2.0,
            w/2.0, w/2.0, w/2.0, w/2.0, w/2.0, w/2.0, w/2.0, w/2.0, w/2.0).finished();

  /**
   * Update step init
   */
  R_lzr_ = MatrixXd::Zero(2, 2);
  R_lzr_(0,0) = std_laspx_*std_laspx_;
  R_lzr_(1,1) = std_laspy_*std_laspy_;
  Zsig_lzr_ = MatrixXd(2, 2 * n_aug_ + 1);
  z_pred_lzr_ = VectorXd(2);
  z_raw_lzr_ = VectorXd(2);
  S_lzr_ = MatrixXd(2, 2);
  Tc_lzr_ = MatrixXd(n_x_, 2);

  R_rdr_ = MatrixXd::Zero(3, 3);
  R_rdr_(0,0) = std_radr_*std_radr_;
  R_rdr_(1,1) = std_radphi_*std_radphi_;
  R_rdr_(2,2) = std_radrd_*std_radrd_;
  Zsig_rdr_ = MatrixXd(3, 2 * n_aug_ + 1);
  z_pred_rdr_ = VectorXd(3);
  z_raw_rdr_ = VectorXd(3);
  S_rdr_ = MatrixXd(3, 3);
  Tc_rdr_ = MatrixXd(n_x_, 3);

  // NIS
  /*NIS_rdr_strm_.open("NIS_rdr_out.txt");
  NIS_lzr_strm_.open("NIS_lzr_out.txt");*/
}

UKF::~UKF() {
    /*NIS_rdr_strm_.close();
    NIS_lzr_strm_.close();*/
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */

  if (!is_initialized_) {
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
          // Convert radar from polar to cartesian coordinates
          //         and initialize state.
          VectorXd cartV = this->tools.ConvertPolarToCartesian(meas_package.raw_measurements_);
          x_ = cartV;
      }
      else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
          // Initialize state.
          x_ << meas_package.raw_measurements_[0],
                  meas_package.raw_measurements_[1],
                  0,
                  0,
                  0;
      }


      previous_timestamp_ = meas_package.timestamp_;
      is_initialized_ = true;
      return;
  }

  double dt = (meas_package.timestamp_ - previous_timestamp_)/ 1000000.0;
  previous_timestamp_ = meas_package.timestamp_;

  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    if (use_radar_)
        UpdateRadar(meas_package);
  } else {
    if (use_laser_)
        UpdateLidar(meas_package);
  }
}

void UKF::Prediction(double dt) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */

  // Generate AugmentedSigmaPoints
  // --------------------------------
    P_aug_ << P_, MatrixXd::Zero(5,2), MatrixXd::Zero(2,5), Q_;

    // create augmented mean vector
    VectorXd x_aug = (VectorXd(7)  << x_, 0, 0).finished();

    MatrixXd tmp = (MatrixXd(7, 7) <<  x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug).finished();
    MatrixXd A = P_aug_.llt().matrixL();

    // set sigma points as columns of matrix Xsig
    Xsig_aug_ << x_aug, tmp + sqrt(lambda_ + n_aug_) * A, tmp - sqrt(lambda_ + n_aug_) * A;

    // SigmaPointPrediction: predict sigma points
    // --------------------------------
    double dt2 = (dt*dt)/2.;
    for (int j=0; j< 2*n_aug_ + 1; ++j) {
        VectorXd col = Xsig_aug_.col(j);
        VectorXd sig = col.head(n_x_);

        double v = col(2),
                psi = col(3),
                psi_d = col(4),
                ua = col(5),
                uf = col(6);

        // noise factor
        VectorXd du = (VectorXd(5) << dt2*cos(psi)*ua,
                dt2*sin(psi)*ua,
                dt*ua,
                dt2 * uf,
                dt * uf).finished();

        VectorXd dx = VectorXd(5);
        if (fabs(psi_d) > 0.001)
            dx << (v/psi_d) * (sin(psi + psi_d*dt) - sin(psi)),
                    (v/psi_d) * (-cos(psi+psi_d*dt) + cos(psi)),
                    0, psi_d * dt, 0;
        else
            dx << v*cos(psi)*dt, v*sin(psi)*dt, 0, psi_d * dt, 0;

        Xsig_pred_.col(j) = sig + dx + du;
    }

    // PredictMeanAndCovariance
    // --------------------------------
    // predict state mean
    x_ = weights_.transpose() * Xsig_pred_.transpose();

    // predict state covariance matrix
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred_.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }
}

/**
 * Update Lidar measurements
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    // PredictLidarMeasurement
    // ---------------------------
    // transform sigma points into measurement space
    Zsig_lzr_.fill(0.0);
    for (int i=0; i < 2 * n_aug_ + 1; ++i) {
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);

        Zsig_lzr_(0,i) = px;
        Zsig_lzr_(1,i) = py;
    }
    // calculate mean predicted measurement
    z_pred_lzr_ =  weights_.transpose() * Zsig_lzr_.transpose();

    // calculate innovation covariance matrix S_lzr_
    S_lzr_.fill(0.0);
    for (int k=0; k < 2*n_aug_+1; ++k) {
        VectorXd z_diff = Zsig_lzr_.col(k) - z_pred_lzr_;
        S_lzr_ += weights_(k) * z_diff * z_diff.transpose();
    }
    // radar measurement noise
    S_lzr_ += R_lzr_;

    z_raw_lzr_ << meas_package.raw_measurements_(0),
            meas_package.raw_measurements_(1);

    // calculate cross correlation matrix
    Tc_lzr_.fill(0.0);
    for (int c=0; c<2*n_aug_+1; ++c) {
        VectorXd x_diff = Xsig_pred_.col(c) - x_;
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        VectorXd z_diff = Zsig_lzr_.col(c) - z_pred_lzr_;

        Tc_lzr_ += weights_(c) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K;
    MatrixXd Sinv = S_lzr_.inverse();
    MatrixXd K = Tc_lzr_ * Sinv;

    // update state mean and covariance matrix
    VectorXd z_diff = z_raw_lzr_ - z_pred_lzr_;
    x_ += K*z_diff;
    P_ -= K*S_lzr_*K.transpose();

    // record NIS values to be plotted later
    /*double nis = z_diff.transpose() * Sinv * z_diff;
    NIS_lzr_strm_ << nis << endl;
    if (nis > 5.991)
        NIS_cnt_gt_chi_lzr_++;
    NIS_cnt_lzr_++;

    cout << "Lidar % above chi_sq (5.991)" << NIS_cnt_gt_chi_lzr_*1.0/NIS_cnt_lzr_ << endl;*/
}

/**
 * Update Radar measurements
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
    // PredictRadarMeasurement
    // ---------------------------
    // transform sigma points into measurement space
    Zsig_rdr_.fill(0.0);
    for (int i=0; i < 2 * n_aug_ + 1; ++i) {
        double px = Xsig_pred_(0,i);
        double py = Xsig_pred_(1,i);
        double v = Xsig_pred_(2,i);
        double psi = Xsig_pred_(3,i);

        double ro = sqrt(px*px + py*py);
        ro = (ro == 0 ? 0.00000001 : ro);
        Zsig_rdr_(0,i) = ro;
        double phi = atan2(py, px);
        Zsig_rdr_(1,i) = phi;
        Zsig_rdr_(2,i) = (px*cos(psi)*v + py*sin(psi)*v)/ro;
    }
    // calculate mean predicted measurement
    z_pred_rdr_ =  weights_.transpose() * Zsig_rdr_.transpose();

    // calculate innovation covariance matrix S
    S_rdr_.fill(0.0);
    for (int k=0; k < 2*n_aug_+1; ++k) {
        VectorXd z_diff = Zsig_rdr_.col(k) - z_pred_rdr_;
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
        S_rdr_ += weights_(k) * z_diff * z_diff.transpose();
    }
    // radar measurement noise
    S_rdr_ += R_rdr_;

    z_raw_rdr_ << meas_package.raw_measurements_(0),
                meas_package.raw_measurements_(1),
                meas_package.raw_measurements_(2);

    // calculate cross correlation matrix
    Tc_rdr_.fill(0.0);
    for (int c=0; c<2*n_aug_+1; ++c) {
        VectorXd x_diff = Xsig_pred_.col(c) - x_;
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
        VectorXd z_diff = Zsig_rdr_.col(c) - z_pred_rdr_;
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        Tc_rdr_ += weights_(c) * x_diff * z_diff.transpose();
    }

    // calculate Kalman gain K;
    MatrixXd Sinv = S_rdr_.inverse();
    MatrixXd K = Tc_rdr_ * Sinv;

    // update state mean and covariance matrix
    VectorXd z_diff = z_raw_rdr_ - z_pred_rdr_;
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    x_ += K*z_diff;
    P_ -= K*S_rdr_*K.transpose();

    // record NIS values to be plotted later
/*    double nis = z_diff.transpose() * Sinv * z_diff;
    NIS_rdr_strm_ << nis << endl;
    if (nis > 7.815)
        NIS_cnt_gt_chi_rdr_++;
    NIS_cnt_rdr_++;

    cout << "Radar % above chi_sq (7.815)" << NIS_cnt_gt_chi_rdr_*1.0/NIS_cnt_rdr_ << endl;*/
}