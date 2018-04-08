#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_  << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

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
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      double rho = measurement_pack.raw_measurements_[0]; // range
      double phi = measurement_pack.raw_measurements_[1]; // bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // velocity of rho

      // convert from polar to cartesian
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = 0;
      double vy = 0;

      ekf_.x_ << px, py, vx, vy;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    // done initializing, no need to predict or update
    if (fabs(ekf_.x_(0)) < 0.0001 and fabs(ekf_.x_(1)) < 0.0001){
      ekf_.x_(0) = 0.0001;
      ekf_.x_(1) = 0.0001;
    }

    // Initial covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

    cout << "EKF init: " << ekf_.x_ << endl;

    // Save the initial timestamp to calc dt later
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  float dt = (measurement_pack.timestamp_ - previous_timestamp_);
  dt /= 1000000.0; // convert micro-second to second
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
          0, 1, 0, dt,
          0, 0, 1, 0,
          0, 0, 0, 1;

  // Noise covariance matrix
  float noise_ax = 9.0;  // per required in task
  float noise_ay = 9.0;  // per required in task

  // Compute covariance matrix Q
  float dt_2 = dt * dt;  // dt^2
  float dt_3 = dt_2 * dt;  // dt^3
  float dt_4 = dt_3 * dt;  // dt^4
  float dt_4_4 = dt_4 / 4;  // dt^4/4
  float dt_3_2 = dt_3 / 2;  // dt^3/2

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_4 * noise_ax, 0, 		dt_3_2 * noise_ax, 	0,
          0, 		dt_4_4 * noise_ay, 0, 			dt_3_2 * noise_ay,
          dt_3_2 * noise_ax, 0, 		dt_2 * noise_ax, 	0,
          0, 		dt_3_2 * noise_ay, 0, 			dt_2 * noise_ay;

  // Check if dt is above a certain threshold before predicting
  if (dt >= 0.000001){
    ekf_.Predict();
  }

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // Jacobian matrix
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

}
