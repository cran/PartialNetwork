#ifndef PACKAGECPdetect_UPD_H
#define PACKAGECPdetect_UPD_H

#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include "vMF.h"


using namespace Rcpp;

void updGnorm (List& Gnorm,
               const List& prior,
               const List& ListIndex,
               const Rcpp::IntegerVector& N,
               const int& M,
               const List& y,
               List& A,
               List& Ay,
               const List& Xb,
               const List& Xgamma,
               const double& alpha,
               const double& sigma2);
void updGnormnoc (List& Gnorm,
                  const List& prior,
                  const List& ListIndex,
                  const Rcpp::IntegerVector& N,
                  const int& M,
                  const List& y,
                  List& A,
                  List& Ay,
                  const List& Xb,
                  const double& alpha,
                  const double& sigma2);
void updGnormblock (List& Gnorm,
                    const List& prior,
                    const List& ListIndex,
                    const Rcpp::IntegerVector& N,
                    const int& M,
                    const List& y,
                    List& A,
                    List& Ay,
                    const List& Xb,
                    const List& Xgamma,
                    const double& alpha,
                    const double& sigma2,
                    const int& nupmax);
void updGnormblocknoc (List& Gnorm,
                       const List& prior,
                       const List& ListIndex,
                       const Rcpp::IntegerVector& N,
                       const int& M,
                       const List& y,
                       List& A,
                       List& Ay,
                       const List& Xb,
                       const double& alpha,
                       const double& sigma2,
                       const int& nupmax);
void updtheta (arma::vec& theta,
               List& Vtheta,
               List& Xb,
               List& Xgamma,
               const double& sigma2,
               const double& kv,
               const double& kbeta,
               const List& Xone,
               const List& X,
               const List& Ay,
               const List& V,
               const arma::vec& invsigmathetatheta0,
               const arma::mat& invsigmatheta,
               const double& M);
void updthetanoc (arma::vec& theta,
                  List& Vtheta,
                  const double& sigma2,
                  const double& kv,
                  const List& Ay,
                  const List& V,
                  const arma::vec& invsigmathetatheta0,
                  const arma::mat& invsigmatheta,
                  const double& M);
void updsigma2 (double& sigma2,
                const arma::vec& theta,
                const double& a,
                const double& b,
                const arma::vec theta0,
                const arma::mat& invsigmatheta,
                const List& Ay,
                const List& Vtheta,
                const double& sumN,
                const double& M);
void updzeta (double& zeta,
              double& alpha,
              List& A,
              double& sumlogdetA,
              List& Ay,
              const List& Gnorm,
              const List& y,
              const double& sigma2,
              const List& Vtheta,
              const double& jumpzeta,
              double& zetaaccept,
              const double& zeta0,
              const double& invsigma2zeta,
              const Rcpp::IntegerVector N,
              const double M);

void updrhopl(List& Gnorm,
              List& prior,
              List& G0obs,
              List& ListIndex,
              arma::vec& rho,
              Eigen::VectorXd& lFdZrhoE1,
              Eigen::VectorXd& lFdZrhoE0,
              const Rcpp::NumericVector weight,
              const arma::mat& dZ,
              const arma::vec& murho,
              const arma::mat& iVrho,
              const arma::mat& jumprho,
              const int& Krho,
              const Rcpp::IntegerVector& N,
              const int& M,
              double& rhoaccept,
              const int& type,
              const bool& Afixed,
              const Eigen::ArrayXd& G0obsvec);

void updrhoARD(List& Gnorm,
               List& prior,
               List& G0obs,
               List& ListIndex,
               List& rho,
               const List& d,
               const arma::vec& zeta,
               const List& murho,
               const List& iVrho,
               const List& jumprho,
               const arma::vec& Krho,
               List& neighbor,
               List& weight,
               List& iARD,
               List& inonARD,
               const Rcpp::IntegerVector& N,
               const Rcpp::IntegerVector& N1,
               const int& M,
               const Rcpp::IntegerVector& P,
               arma::vec& rhoaccept,
               const arma::vec& type);

#endif