#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include "tools.h"
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:
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
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // noise covariance matrix
  Eigen::MatrixXd Q_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  const int n_x_ = 5;

  // Augmented state dimension
  const int n_aug_ = 7;

  // Sigma point spreading parameter
  const double lambda_ = 3 - n_x_;

  long long int previous_timestamp_ = 0;

  Tools tools;

  // Prediction
  MatrixXd Xsig_aug_;  // Augmented Sigma points matrix
  MatrixXd P_aug_;     // Augmented Covariance matrix

  //Measurement Update
  Eigen::MatrixXd R_rdr_;
  MatrixXd Zsig_rdr_;   // matrix for sigma points in measurement space
  VectorXd z_pred_rdr_;  // mean predicted measurement
  VectorXd z_raw_rdr_;   //raw measurements
  MatrixXd S_rdr_;    // measurement covariance matrix S
  MatrixXd Tc_rdr_; //matrix for cross correlation Tc

  Eigen::MatrixXd R_lzr_;
  MatrixXd Zsig_lzr_;
  VectorXd z_pred_lzr_;
  VectorXd z_raw_lzr_;
  MatrixXd S_lzr_;
  MatrixXd Tc_lzr_;

  // NIS analysis
  std::ofstream NIS_rdr_strm_;
  int NIS_cnt_rdr_ = 0;
  int NIS_cnt_gt_chi_rdr_ = 0;

  std::ofstream NIS_lzr_strm_;
  int NIS_cnt_lzr_ = 0;
  int NIS_cnt_gt_chi_lzr_ = 0;

/*    // Analyze possible max values for tracked object acceleration and yaw acceleration
  double prev_v_ = 0.;
  double prev_yawd_ = 0.;
  double acc_a = 0.;
  double acc_ydd = 0.;
  double acc_a_max = 0.;
  double acc_ydd_max = 0.;*/
};

#endif  // UKF_H