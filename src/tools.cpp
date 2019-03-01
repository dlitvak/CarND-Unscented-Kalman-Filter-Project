#include "tools.h"

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

Eigen::VectorXd Tools::ConvertPolarToCartesian(const Eigen::VectorXd &polarV) {
    VectorXd cart = VectorXd(5);
    double ro = polarV(0);
    double phi = polarV(1);

    double px = ro * cos(phi);
    double py = ro * sin(phi);

    // We can't determine cartesian velocity from roDot.  Setting it to 0
    cart << px, py, 0, 0, 0;
    return cart;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
    VectorXd rmse(4);
    rmse << 0,0,0,0;

    size_t esz = estimations.size();
    size_t gsz = ground_truth.size();
    if (esz == 0 || esz != gsz) {
        return rmse;
    }

    for (unsigned int i=0; i < estimations.size(); ++i) {
        VectorXd diffV = estimations[i] - ground_truth[i];
        diffV = diffV.array()*diffV.array();

        rmse += diffV;
    }

    rmse /= estimations.size();
    rmse = rmse.array().sqrt();
    return rmse;
}