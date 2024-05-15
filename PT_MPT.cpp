#include <Rcpp.h>
#include <cmath>
#include <math.h>
using namespace Rcpp;

// Need helpder function for product of a vector
double prod_vec(NumericVector vec) {
  double ans = 1.0;
  int n = static_cast<int>(vec.length());
  for (int i = 0; i < n; i++) {
    ans *= vec[i];
  }
  return ans;
}

// We first need a few helper function to help working with PT and MPT
// This is in a function to help finding the range of a datapoint in the number line (for PT based on Normal distribution)
// [[Rcpp::export]]
int find_quantile_cpp(double x, double mu, double sigma, int n_quantiles) {
  // need a variable to store the answer
  int ans = -1;
  // Conduct binary search on this
  int l = 1; 
  int r = n_quantiles;
  int m;
  while (l <= r) {
    m = (l + r) / 2;
    if (R::qnorm((double)m / n_quantiles, mu, sigma, true, false) >= x) {
      ans = m;
      r = m - 1;
    } else {
      l = m + 1;
    }
  }
  return ans;
}

// We also need another function to generate theta
// We store it as a matrix => less need of creating new items
// [[Rcpp::export]]
NumericMatrix draw_theta_cpp(bool has_prev, IntegerMatrix data_count, int c, int J) {
  // data_count and theta matrix both have: J * 2^J
  // store the resultng theta
  NumericMatrix sampled_theta(J, (int)std::pow(2, J));
  // variable for iteration
  int k,l;
  // now we get each theta
  for (k = 0; k < J; k++) {
    // level k => 2^(k+1) value
    // odd value from 0 -> 2^(k+1) - 2 <=> 0 -> 2^(k) - 1
    for (l = 0; l < (int)std::pow(2, k); l++) {
      // make the new theta
      // add to the position in resulting matrix
      if (has_prev) {
        // consider the previous data since we are working with posterior dist now
        sampled_theta(k, 2*l) = R::rbeta(data_count(k,2*l) + c*std::pow(k, 2), data_count(k,2*l + 1) + c*std::pow(k, 2));
      } else {
        sampled_theta(k, 2*l) = R::rbeta(c*std::pow(k, 2), c*std::pow(k, 2));
      }
      sampled_theta(k, 2*l + 1) = 1 - sampled_theta(k, 2*l);
    }
    
  }
  return sampled_theta;
}

// IMPLEMENTATION OF A RANDOM VARIABLE DRAWN FROM A POLYA TREE PRIOR
// We first need the pdf for Polya tree given the normal distribution and the thetas
// [[Rcpp::export]]
double pdf_pt_prior_cpp(double x, double mu, double sigma, NumericMatrix theta_list) {
  // mu: mean of the normal distribution
  // sigma: sd of the normal distribution
  // theta_list: a list of thetas for J level - presented as a matrix of J * 2^J, need to consider how many values
  //             are the actual thetas when working with this
  // We first construct the level of Polya tree, J
  int J = theta_list.nrow();
  // Now we need to check the quantile of x (if we divide the normal distribution to 2^J quantile)
  int curr_quantile = find_quantile_cpp(x, mu, sigma, (int)std::pow(2, J));
  // After we found the quantile of x, we need to find the theta that associated with it so that we can multiply
  // We do this by keep dividing quantile_x by 2, until we get J value
  double theta_prod = 1.0;
  for (int i = J - 1; i >= 0; i--) {
    // We get the current quantile 
    theta_prod *= theta_list(i, curr_quantile - 1); // need to -1 since we are working with 0-index list
    curr_quantile = (curr_quantile / 2) + (curr_quantile % 2);
  }
  return R::dnorm(x, mu, sigma, false) * std::pow(2, J) * theta_prod;
}

// Now we make cdf based on the pdf, we will do this based on simple calculation of riemann sum with many terms
// [[Rcpp::export]]
double cdf_pt_prior_cpp(double x, double mu, double sigma, NumericMatrix theta_list) {
  // mu: mean of the normal distribution
  // sigma: sd of the normal distribution
  // theta_list: a list of thetas for J level - presented as a matrix of J * 2^J, need to consider how many values
  //             are the actual thetas when working with this
  // in order to approximate using riemann sum, we can consider the range of [-5sd, x] (use -5sd instead of -Inf)
  // we iter and calculate the mid value directly
  double inc = (x - mu + 5*sigma) / 100000;
  // store current mid value for the integration
  double curr_x = mu - 5*sigma + inc / 2;
  // store the answer
  double ans = 0.0;
  for (int i = 0; i < 100000; i++) {
    ans += pdf_pt_prior_cpp(curr_x, mu, sigma, theta_list);
    curr_x += inc;
  }
  return ans * inc;
}

// Implementation of pdf and cdf of a mixture distribution over theta with Polya tree prior
// (since theta will be drawn randomly based on c, so we will take a mixture over theta)
// We first need to implement an pdf for this case of polya tree (with more random variable)
// X | c, mu, sigma = integral of X | theta, c, mu, sigma * theta | c dtheta
// Since there can be many value of theta, we can consider sampling the theta
// [[Rcpp::export]]
double pdf_pt_prior_mixture_c_cpp(double x, double mu, double sigma, 
                                  int c, int J, int n_iter, 
                                  bool has_prev, IntegerMatrix data_count) {
  // mu, sigma: mean and sd of the normal distribution
  // c: constant in the distribution for theta
  // J: number of level
  // data_count: count of data at each quantile at each level of tree 
  // has_orev: whether our distribution has had past data
  // Firstly, noticing that all of our theta are in the range [0, 1], so the volume space of all theta is finite, so we can use a Monte Carlo integration method
  // Now we will sample theta based on INDEPENDENT distribution of the thetas
  // and for each set, we will calculate the pdf (using the function that we madee
  double ans = 0.0;
  // need a variable for the sampled theta
  NumericMatrix sampled_theta;
  for (int i = 0; i < n_iter; i++) {
    // sample our thetas and save it in a list
    sampled_theta = draw_theta_cpp(has_prev, data_count, c, J);
    // Now we calculate integral of X | theta, c, mu, sigma * p(theta) = E(X | theta, c, mu, sigma) => 
    // we calculate the mean by each turn calculate the pdf
    // and then we add the result to ans
    ans += pdf_pt_prior_cpp(x, mu, sigma, sampled_theta);
  }
  return ans / n_iter;
}

// We can also make a cdf function based on this
// [[Rcpp::export]]
double cdf_pt_prior_mixture_c_cpp(double x, double mu, double sigma, 
                                  int c, int J, int n_iter, 
                                  bool has_prev, IntegerMatrix data_count) {
  // we iter and calculate the mid value directly
  double inc = (x - mu + 5*sigma) / 10000;
  // store current mid value for the integration
  double curr_x = mu - 5*sigma + inc / 2;
  // store the answer
  double ans = 0.0;
  for (int i = 0; i < 10000; i++) {
    ans += pdf_pt_prior_mixture_c_cpp(curr_x, mu, sigma, c, J, 10000, has_prev, data_count);
    curr_x += inc;
  }
  return ans * inc;
}

// Given the value of c (for distribution of theta) to begin with, we can update posterior by just finding the list of data count for each quantile, and then use the above function for pdf
// [[Rcpp::export]]
IntegerMatrix pt_update_cpp(NumericVector x_vec, double mu_0, double sigma_0, int J) {
  // x_vec: vector of new data
  // mu_0: initial value for mean
  // sigma_0: initial value for sd
  // J: tree level
  // First we need a tree level
  // first of all, for easier update of theta, we need to count how many values in each quantiles, for different number of tree level
  IntegerMatrix quantile_data_count(J, (int)std::pow(2, J)); 
  // store the length of x
  int length_x = static_cast<int>(x_vec.length());
  // store the quantile found
  int curr_quantile;
  // variable for iteration
  int i, j;
  // We do this backwardly from 2^J quantiles to 2 quantiles, since we can use the larger list to create smaller lsit
  for (i = J - 1; i >= 0; i--) {
    if (i == J - 1) {
      // Now we get the position of each data point in 1 of 2^J quantiles
      // and then update it to the count
      for (j = 0; j < length_x; j++) {
        curr_quantile = find_quantile_cpp(x_vec[j], mu_0, sigma_0, (int)std::pow(2, i + 1));
        quantile_data_count(i, curr_quantile - 1) += 1;
      }
    } else {
      for (j = 0; j < (int)std::pow(2, i+1); j++) {
        // formula for updating 
        quantile_data_count(i, j) = quantile_data_count(i + 1, 2*j) + quantile_data_count(i + 1, 2*j + 1);
      }
    }
  }
  // After findind the number of data in each quantiles for each level, we will now update the data to use for pdf
  return quantile_data_count;
}

// Implementation of random variable drawn from a Mixture of Polya Tree 
// (assuming list of theta has already been drawn)
// [[Rcpp::export]]
double pdf_mpt_prior_cpp(double x, 
                         double mu_0, Function d_mu, Function r_mu, // initial mu and distribution of mu
                         double sigma_0, Function d_sigma, Function r_sigma, // initial theta and distribution of theta
                         NumericMatrix theta_list, // list of theta that will be drawn 
                         int n_iter, // number of mu, theta to be drawn
                         bool has_prev, NumericVector prev_data,
                         double t_mu = 1, double t_sigma = 1) {
  // First we need to save the answer
  double ans = 0;
  // and we need to preprare the initial value for mu, sigma
  double prev_mu = mu_0;
  double prev_sigma = sigma_0;
  // need variable to store candidate value
  double candidate_mu;
  double candidate_sigma;
  // store the rate of acceptance
  double a;
  // store the simulated probability
  double prob_drawn;
  // get the length of previous data
  int n = static_cast<int>(prev_data.length());
  // For each turn, we draw the mu and sigma based on Metropolis-Hastings (for posterior), or based on their prior 
  for (int i = 0; i < n_iter; i++) {
    if (has_prev) {
      // Draw based on Metropolis-Hasting since this is a posterior distribution
      // We first draw mu
      candidate_mu = R::rnorm(prev_mu, prev_sigma * sqrt(t_mu / n));
      // calculate a to find acceptance criteria 
      a = prod_vec(dnorm(prev_data, candidate_mu, prev_sigma)) * as<double>(d_mu(candidate_mu)) / (prod_vec(dnorm(prev_data, prev_mu, prev_sigma)) * as<double>(d_mu(prev_mu)));
      if (a <= 1) {
        prob_drawn = R::runif(0, 1);
        if (prob_drawn > a) {
          candidate_mu = prev_mu; 
        }
      }
      // We then draw sigma
      candidate_sigma = sqrt(R::rlnorm(log(std::pow(prev_sigma, 2)), sqrt(t_sigma)));
      // calculate a for parameter update
      a = prod_vec(dnorm(prev_data, prev_mu, candidate_sigma)) * as<double>(d_sigma(candidate_sigma)) / (prod_vec(dnorm(prev_data, prev_mu, prev_sigma)) * as<double>(d_sigma(prev_sigma)));
      if (a <= 1)  {
        prob_drawn = R::runif(0, 1);
        if (prob_drawn > a) {
          candidate_sigma = prev_sigma;
        }
      }
    } else {
      candidate_mu = as<double>(r_mu());
      candidate_sigma = as<double>(r_sigma());
    }
    ans = ans + pdf_pt_prior_cpp(x, candidate_mu, candidate_sigma, theta_list);
  }
  return ans/n_iter;
}
