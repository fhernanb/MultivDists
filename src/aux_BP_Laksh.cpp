#include <Rcpp.h>
using namespace Rcpp;

#include <Rcpp.h>
#include <cmath>

//' Auxiliar cpp function to obtain the probability for a vector x.
//' @param input numeric value with the information.
//' @export
//' @return returns the pmf for a vector.
// [[Rcpp::export]]
NumericVector aux_BP_Laksh(NumericVector input) {
  double x1 = input[0];
  double x2 = input[1];
  double l1 = input[2];
  double l2 = input[3];
  double alpha = input[4];

  double k = 1 - std::exp(-1);
  double A = std::exp(-l1 * k);
  double B = std::exp(-l2 * k);

  double part1 = (-l1 - l2) + x1 * std::log(l1) + x2 * std::log(l2) - std::lgamma(x1 + 1) - std::lgamma(x2 + 1);
  double part2 = std::log(1 + alpha * (std::exp(-x1) - A) * (std::exp(-x2) - B));

  return NumericVector::create(part1 + part2);
}

