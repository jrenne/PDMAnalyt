
// ---- Enable C++11: ----
// [[Rcpp::plugins("cpp11")]]

#include <Rcpp.h>
#include <RcppEigen.h>
#include <unsupported/Eigen/SpecialFunctions>

//[[Rcpp::depends(RcppEigen)]] 

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
Rcpp::List f_essai_cpp(const  int & i,
                       const  double & x){
  
  double i2 = pow(i,2) ;
  
  return List::create(Named("i2")  = i2) ;
  
}



Eigen::MatrixXd na_matrix(int r,
                          int c){
  // create matrix with NaNs
  Eigen::MatrixXd m(r,c) ;
  m = MatrixXd::Zero(r,c) ;
  //m(0,0) = R_NaN ;
  m.fill(R_NaN) ;
  return m ;
}

// Eigen::VectorXd getEigenValues(const Eigen::MatrixXd & M) {
//   SelfAdjointEigenSolver<MatrixXd> es(M);
//   return es.eigenvalues();
// }
// 
// // [[Rcpp::export]]
// Rcpp::List Eigenvalue_decomp(const MatrixXd A){
//   
//   Eigen::EigenSolver<MatrixXd> eig(A);
//   Eigen::MatrixXcd Vec = eig.eigenvectors();
//   Eigen::VectorXcd Val = eig.eigenvalues();
//   
//   return(
//     List::create(Named("values")  = Val,
//                  Named("vectors") = Vec)
//   );
//   
// }

Eigen::MatrixXd which_min_row(const Eigen::MatrixXd & M){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd VecMins = Eigen::MatrixXd::Zero(p,1) ;
  for(int i = 0; i < p; i++){
    int aux = 0 ;
    double mini = M(i,0) ;
    double v = 0 ;
    for(int j = 1; j < n; j++){
      v = M(i,j) ;
      if( v < mini ){
        aux = j ;
        mini = v ;
      }
    }
    VecMins(i) = aux + 1;
  }
  return VecMins;
}

Eigen::MatrixXd find_closest_in_vec(const Eigen::MatrixXd & M,
                                    const Eigen::MatrixXd & vec){
  int p = M.rows() ;
  int n = M.cols() ;
  int vec_values = vec.rows() ;
  Eigen::MatrixXd MatIndicators = Eigen::MatrixXd::Zero(p,n) ;
  double smallest_dist ;
  double dist ;
  double coeff ;
  double coeff_vec ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      coeff = M(i,j) ;
      coeff_vec = vec(0,0) ;
      smallest_dist = (coeff - coeff_vec) * (coeff - coeff_vec) ;
      MatIndicators(i,j) = 1 ; // may stay there if no one smaller
      for(int k = 1; k < vec_values; k++){
        coeff_vec = vec(k,0) ;
        dist = (coeff - coeff_vec) * (coeff - coeff_vec) ;
        if(dist < smallest_dist){
          MatIndicators(i,j) = k + 1;
          smallest_dist = dist ;
        }
      }
    }
  }
  return MatIndicators ;
}

Eigen::MatrixXd fill_from_indic(const Eigen::MatrixXd & MatIndic,
                                const Eigen::MatrixXd & vec){
  int p = MatIndic.rows() ;
  int n = MatIndic.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,n) ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      int k = MatIndic(i,j) - 1;
      MatOutput(i,j) = vec(k,0) ;
    }
  }
  return MatOutput;
}



Eigen::MatrixXd mult(const Eigen::MatrixXd & matrix,
                     const double & scalar){
  int p = matrix.rows() ;
  int n = matrix.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,n) ;
  Eigen::MatrixXd AUX = Eigen::MatrixXd::Constant(p, n, scalar) ;
  MatOutput = matrix.array() * AUX.array() ;
  return MatOutput;
}

Eigen::MatrixXd add(const Eigen::MatrixXd & matrix,
                    const double & scalar){
  int p = matrix.rows() ;
  int n = matrix.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,n) ;
  Eigen::MatrixXd AUX = Eigen::MatrixXd::Constant(p, n, scalar) ;
  MatOutput = matrix + AUX ;
  return MatOutput;
}


// [[Rcpp::export]]
Eigen::MatrixXd mean_per_row(const Eigen::MatrixXd & M){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd MatOutput = Eigen::MatrixXd::Zero(p,1) ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      MatOutput(i,0) = MatOutput(i,0) + 1/(n*1.0) * M(i,j) ;
    }
  }
  return MatOutput;
}

Eigen::MatrixXd pmax_cpp(const Eigen::MatrixXd & M,
                         double x){
  int p = M.rows() ;
  int n = M.cols() ;
  Eigen::MatrixXd MatOutput = M ;
  for(int i = 0; i < p; i++){
    for(int j = 0; j < n; j++){
      if(M(i,j)<x){
        MatOutput(i,j) = x ;
      }
    }
  }
  return MatOutput;
}


Eigen::MatrixXd seqEigen(int start,
                         int end){
  int length = end - start + 1 ;
  Eigen::MatrixXd V = Eigen::MatrixXd::Zero(length,1) ;
  for(int i = 0; i < length; i++){
    V(i,0) = i + start ;
  }
  return V ;
}

// [[Rcpp::export]]
Eigen::MatrixXd Mpower(const Eigen::MatrixXd & A,
                       int   n){
  // compute A^(2^n)
  int rA = A.rows() ;
  int cA = A.cols() ;
  
  Eigen::MatrixXd M   = Eigen::MatrixXd::Zero(rA, cA) ;
  M = A ;
  
  for(int i = 1; i <= n; i++){
    M = M * M ;
  }
  return M ;
}

Eigen::MatrixXd kronecker_cpp(const Eigen::MatrixXd & A,
                              const Eigen::MatrixXd & B){
  int rA = A.rows() ;
  int cA = A.cols() ;
  int rB = B.rows() ;
  int cB = B.cols() ;
  
  Eigen::MatrixXd M   = Eigen::MatrixXd::Zero(rA*rB, cA*cB) ;
  Eigen::MatrixXd AUX = Eigen::MatrixXd::Zero(rB, cB) ;
  
  for(int i = 1; i <= rA; i++){
    for(int j = 1; j <= cA; j++){
      AUX = Eigen::MatrixXd::Constant(rB, cB, A(i-1,j-1)) ;
      AUX = AUX.array() * B.array() ;
      M.block((i-1)*rB,(j-1)*cB,rB,cB) = AUX ;
    }
  }
  return M ;
}

Eigen::MatrixXd Vec(const Eigen::MatrixXd & M){
  int rM = M.rows() ;
  int cM = M.cols() ;
  Eigen::MatrixXd vecM  = Eigen::MatrixXd::Zero(rM*cM, 1);
  for(int j = 1; j <= cM; j++){
    vecM.block((j-1)*rM,0,rM,1) = M.col(j-1) ;
  }
  return vecM ;
}

Eigen::MatrixXd Vec_1(const Eigen::MatrixXd & vecM){
  int n2 = vecM.rows() ;
  int n = sqrt(n2) ;
  
  Eigen::MatrixXd M  = Eigen::MatrixXd::Zero(n, n);
  for(int j = 1; j <= n; j++){
    M.block(0,j-1,n,1) = vecM.block((j-1)*n,0,n,1) ;
  }
  return M ;
}

Eigen::MatrixXd pnorm_cpp(Eigen::MatrixXd X){
  int nb_r = X.rows() ;
  int nb_c = X.cols() ;
  Eigen::MatrixXd p    = Eigen::MatrixXd::Zero(nb_r, nb_c);
  Eigen::MatrixXd Cst1 = Eigen::MatrixXd::Constant(nb_r, nb_c, 0.5) ;
  Eigen::MatrixXd Cst2 = Eigen::MatrixXd::Constant(nb_r, nb_c, -sqrt(0.5)) ;
  p = Cst1.array() * ((Cst2.array() * X.array()).erfc()).array() ;
  return p;
}

Eigen::MatrixXd dnorm_cpp(Eigen::MatrixXd X){
  double pi = 3.141592653589793 ;
  int nb_r = X.rows() ;
  int nb_c = X.cols() ;
  Eigen::MatrixXd p  = Eigen::MatrixXd::Zero(nb_r, nb_c);
  Eigen::MatrixXd Cst1 = Eigen::MatrixXd::Constant(nb_r, nb_c, - 0.5) ;
  Eigen::MatrixXd Cst2 = Eigen::MatrixXd::Constant(nb_r, nb_c, 1/(sqrt(2) * sqrt(pi))) ;
  p = Cst2.array() * (Cst1.array() * (X.array().square()).array()).exp() ;
  return p;
}

// [[Rcpp::export]]
Eigen::MatrixXd apply_cumsum(Eigen::MatrixXd M){
  int nb_c = M.cols() ;
  Eigen::MatrixXd MatriX ;
  MatriX = M ;
  for (int i = 1; i <= nb_c - 1; ++i){
    MatriX.col(i) = MatriX.col(i-1) + M.col(i) ;
  }
  return MatriX ;
}


// // [[Rcpp::export]]
// Rcpp::List solve_ToyModel_OLD(const Eigen::MatrixXd all_d,
//                               const Eigen::MatrixXd all_rr,
//                               const Eigen::MatrixXd all_eps,
//                               const Rcpp::List Model,
//                               const int nb_iter
// ){
//   
//   double gamma     = Model("gamma") ;
//   double mu        = Model("mu") ;
//   double delta     = Model("delta") ;
//   double b         = Model("b") ;
//   double chi       = Model("chi") ;
//   double beta      = Model("beta") ;
//   double d_star    = Model("d_star") ;
//   double sigma_eps = Model("sigma_eps") ;
//   double alpha     = Model("alpha") ;
//   double s_star    = Model("s_star") ;
//   double RR        = Model("RR") ;
//   
//   int nb_grid_d  = all_d.rows() ;
//   int nb_grid_rr = all_rr.rows() ;
//   int nb_states  = nb_grid_d * nb_grid_d * nb_grid_rr ;
//   int nb_eps     = all_eps.rows() ;
//   
//   Eigen::MatrixXd vec_1_d   = Eigen::MatrixXd::Constant(nb_grid_d, 1, 1) ;
//   Eigen::MatrixXd vec_1_rr  = Eigen::MatrixXd::Constant(nb_grid_rr, 1, 1) ;
//   Eigen::MatrixXd vec_1_eps = Eigen::MatrixXd::Constant(nb_eps, 1, 1) ;
//   Eigen::MatrixXd vec_1_x   = Eigen::MatrixXd::Constant(nb_states, 1, 1) ;
//   Eigen::MatrixXd mat_1     = Eigen::MatrixXd::Constant(nb_states, nb_eps, 1) ;
//   Eigen::MatrixXd ALL_eps   = Eigen::MatrixXd::Constant(nb_states, nb_eps, 1) ;
//   ALL_eps = vec_1_x * all_eps.transpose() ;
//   
//   Eigen::MatrixXd Mat_1 = Eigen::MatrixXd::Constant(nb_states, nb_eps, 1) ;
//   Eigen::MatrixXd Mat_nb_grid_d = Eigen::MatrixXd::Constant(nb_states, nb_eps, nb_grid_d) ;
//   Eigen::MatrixXd Mat_nb_grid_d2 = Eigen::MatrixXd::Constant(nb_states, nb_eps, nb_grid_d*nb_grid_d) ;
//   
//   Eigen::MatrixXd d   = Eigen::MatrixXd::Zero(nb_grid_d^2*nb_grid_rr, 1) ;
//   Eigen::MatrixXd d_1 = Eigen::MatrixXd::Zero(nb_grid_d^2*nb_grid_rr, 1) ;
//   Eigen::MatrixXd rr  = Eigen::MatrixXd::Zero(nb_grid_d^2*nb_grid_rr, 1) ;
//   d   = kronecker_cpp(kronecker_cpp(vec_1_rr,vec_1_d),all_d) ;
//   d_1 = kronecker_cpp(kronecker_cpp(vec_1_rr,all_d),vec_1_d) ;
//   rr  = kronecker_cpp(kronecker_cpp(all_rr,vec_1_d),vec_1_d) ;
//   
//   double r = - log(delta) + mu ;
//   
//   Eigen::MatrixXd q     = Eigen::MatrixXd::Constant(nb_states, 1, r + .01) ; // vector of sovereign yields
//   Eigen::MatrixXd q0    = Eigen::MatrixXd::Constant(nb_states, 1, r) ; // risk-free yields
//   
//   Eigen::MatrixXd q_iter_1 = q ; // used to check convergence over iterations
//   Eigen::MatrixXd q_chge   = Eigen::MatrixXd::Zero(nb_states, 1) ; // used to check convergence
//   
//   Eigen::MatrixXd rr_tp1    = Eigen::MatrixXd::Zero(nb_states, 1) ;
//   Eigen::MatrixXd d_tp1_aux = Eigen::MatrixXd::Zero(nb_states, 1) ;
//   
//   Eigen::MatrixXd aux_nb_states = Eigen::MatrixXd::Zero(nb_states, 1) ;
//   Eigen::MatrixXd AUX_nb_states = Eigen::MatrixXd::Zero(nb_states, nb_states) ;
//   
//   Eigen::MatrixXd all_d_tp1    = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_d_1_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_rr_tp1   = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   
//   Eigen::MatrixXd indicators_d_tp1    = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd indicators_d_t  = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd indicators_rr_tp1   = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd indicators_x        = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_lambdas         = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_proba_def       = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_q_tp1           = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd all_q0_tp1          = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd E                   = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   Eigen::MatrixXd Q                   = Eigen::MatrixXd::Zero(nb_states, nb_eps) ;
//   
//   Eigen::MatrixXd C1 = Eigen::MatrixXd::Constant(nb_states, 1, chi/(1 + mu)) ;
//   Eigen::MatrixXd C2 = Eigen::MatrixXd::Constant(nb_states, 1, 1/(1 + mu)) ;
//   Eigen::MatrixXd C_beta   = Eigen::MatrixXd::Constant(nb_states, 1, beta) ;
//   Eigen::MatrixXd C_d_star = Eigen::MatrixXd::Constant(nb_states, 1, d_star) ;
//   Eigen::MatrixXd C_s_star = Eigen::MatrixXd::Constant(nb_states, 1, s_star) ;
//   Eigen::MatrixXd C_alpha  = Eigen::MatrixXd::Constant(nb_states, 1, alpha) ;
//   Eigen::MatrixXd C_chi    = Eigen::MatrixXd::Constant(nb_states, nb_eps, chi) ;
//   Eigen::MatrixXd C_xxx    = 
//     Eigen::MatrixXd::Constant(nb_states, nb_eps,
//                               exp(gamma*b)*RR/(1 - chi * delta * exp(- mu))) ;
//   Eigen::MatrixXd C_xxx_RR1    = 
//     Eigen::MatrixXd::Constant(nb_states, nb_eps,
//                               exp(gamma*b)/(1 - chi * delta * exp(- mu))) ;
//   Eigen::MatrixXd C_xxxx   = 
//     Eigen::MatrixXd::Constant(nb_states, nb_eps, exp(-log(delta) + mu)) ;
//   
//   
//   for (int i = 0; i < nb_iter; i++){
//     rr_tp1 = q.array() * (d.array() - C1.array() * d_1.array()) + 
//       C1.array() * rr.array() ;
//     d_tp1_aux  = C2.array() * d.array() - C_beta.array() * (d.array() - 
//       C_d_star.array()) + C2.array() * rr.array() ;
//     
//     all_d_tp1   = d_tp1_aux * vec_1_eps.transpose() - ALL_eps ;
//     all_d_1_tp1 = d * vec_1_eps.transpose() ;
//     all_rr_tp1  = rr_tp1 * vec_1_eps.transpose() ;
//     
//     indicators_d_tp1    = find_closest_in_vec(all_d_tp1,all_d) ;
//     indicators_rr_tp1   = find_closest_in_vec(all_rr_tp1,all_rr) ;
//     
//     indicators_x = indicators_d_tp1.array() +
//       (indicators_d_t - Mat_1).array() * Mat_nb_grid_d.array() +
//       (indicators_rr_tp1 - Mat_1).array() * Mat_nb_grid_d2.array() ;
//     
//     // Compute probabilities of default:
//     aux_nb_states = C_beta.array() * (d - C_d_star).array() ;
//     AUX_nb_states = aux_nb_states * vec_1_eps.transpose() + 
//       ALL_eps - C_s_star * vec_1_eps.transpose();
//     all_lambdas = pmax_cpp(AUX_nb_states, 0) ;
//     all_lambdas = - (C_alpha  * vec_1_eps.transpose()).array() *
//       all_lambdas.array() ;
//     
//     all_proba_def = mat_1.array() - all_lambdas.array().exp() ;
//     
//     // Update q (sovereign yields):
//     all_q_tp1 = fill_from_indic(indicators_x, q) ;
//     AUX_nb_states = (mat_1 + all_q_tp1).array()/(mat_1 + all_q_tp1 - C_chi).array() ;
//     E = AUX_nb_states.array() + all_proba_def.array() * (C_xxx - AUX_nb_states).array() ;
//     Q = C_chi.array() - mat_1.array() + C_xxxx.array() * (E.cwiseInverse()).array() ;
//     q = mean_per_row(Q) ;
//     
//     // Update q0 (risk-free yield, obtained for RR=1):
//     all_q0_tp1 = fill_from_indic(indicators_x, q0) ;
//     AUX_nb_states = (mat_1 + all_q0_tp1).array()/(mat_1 + all_q0_tp1 - C_chi).array() ;
//     E = AUX_nb_states.array() + all_proba_def.array() * (C_xxx_RR1 - AUX_nb_states).array() ;
//     Q = C_chi.array() - mat_1.array() + C_xxxx.array() * (E.cwiseInverse()).array() ;
//     q0 = mean_per_row(Q) ;
//     
//     q_chge = q - q_iter_1 ;
//     q_iter_1 = q; // for next iteration
//   }
//   
//   return List::create(Named("q")  = q,
//                       Named("q0") = q0,
//                       Named("d") = d,
//                       Named("d_1") = d_1,
//                       Named("rr") = rr,
//                       Named("q_chge") = q_chge,
//                       Named("indicators_x") = indicators_x,
//                       Named("all_proba_def") = all_proba_def
//   ) ;
// }


// [[Rcpp::export]]
Rcpp::List compute_SDF(const Rcpp::List Model){
  
  Eigen::MatrixXd Omega  = Model("Omega") ;
  Eigen::MatrixXd mu_pi  = Model("mu_pi") ;
  Eigen::MatrixXd mu_y   = Model("mu_y") ;
  
  double kappa_pi  = Model("kappa_pi") ;
  double kappa_y   = Model("kappa_y") ;
  
  double gamma     = Model("gamma") ;
  double delta     = Model("delta") ;
  
  int nb_m       = Omega.rows() ;
  
  Eigen::MatrixXd vec_1_m       = Eigen::MatrixXd::Constant(nb_m,1,1) ;
  Eigen::MatrixXd mu_u          = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd mu_f0         = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd mu_f1         = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd mu_f1_real    = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd mu_f1_nominal = Eigen::MatrixXd::Zero(nb_m,1) ;
  
  Eigen::MatrixXd aux = Eigen::MatrixXd::Zero(nb_m,1) ;
  
  double delta_gamma = delta/(1 - gamma) ;
  double log_delta   = log(delta) ;
  
  for (int j = 0; j < 1000; j++){
    aux = (mult(mu_u + mu_y,1 - gamma).array()).exp() ;
    aux = aux.transpose() * Omega.transpose() ;
    aux = (aux.array()).log() ;
    mu_u = mult(aux,delta_gamma).transpose() ;
  }
  aux = (mult(mu_u + mu_y,1 - gamma).array()).exp() ;
  aux = aux.transpose() * Omega.transpose() ;
  aux = (aux.array()).log() ;
  mu_f0         = add(-aux,log_delta).transpose() ;
  mu_f1_real    = mult(mu_u,1 - gamma) - mult(mu_y,gamma) ;
  mu_f1_nominal = mu_f1_real - mu_pi ;
  mu_f1         = mu_f1_nominal + mult(mu_pi,kappa_pi) + mult(mu_y,kappa_y) ;
  
  return List::create(Named("mu_u")          = mu_u,
                      Named("mu_f0")         = mu_f0,
                      Named("mu_f1")         = mu_f1,
                      Named("mu_f1_real")    = mu_f1_real,
                      Named("mu_f1_nominal") = mu_f1_nominal) ;
}



// [[Rcpp::export]]
Rcpp::List compute_LTRF_bond_prices(const Rcpp::List Model,
                                    const int maxH){
  
  Eigen::MatrixXd Omega  = Model("Omega") ;
  int nb_m               = Omega.rows() ;
  
  double kappa_pi = Model("kappa_pi") ;
  double kappa_y  = Model("kappa_y") ;
  
  Eigen::MatrixXd mu_pi  = Model("mu_pi") ;
  Eigen::MatrixXd mu_y   = Model("mu_y") ;
  
  Eigen::MatrixXd vec_1_m = Eigen::MatrixXd::Constant(nb_m,1,1) ;
  
  // Get SDF specification:
  Rcpp::List res_SDF = compute_SDF(Model) ;
  Eigen::MatrixXd mu_f0         = res_SDF("mu_f0") ;
  Eigen::MatrixXd mu_f1         = res_SDF("mu_f1") ;
  // Eigen::MatrixXd mu_f1_real    = res_SDF("mu_f1_real") ;
  // Eigen::MatrixXd mu_f1_nominal = res_SDF("mu_f1_nominal") ;
  
  // mu_f1_nominal = mu_f1_nominal.array() + mult(mu_pi,kappa_pi).array() + mult(mu_y,kappa_y).array() ;
  // mu_f1_nominal = mu_f1 ;
  
  Eigen::MatrixXd exp_mu_f1    = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd exp_mu_f1_f0 = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd exp_mu_f0            = Eigen::MatrixXd::Zero(nb_m,1) ;
  
  exp_mu_f1    = (mu_f1.array()).exp() ;
  exp_mu_f1_f0 = ((mu_f1 + mu_f0).array()).exp() ;
  exp_mu_f0            = (mu_f0.array()).exp() ;
  
  Eigen::MatrixXd Mlast = (exp_mu_f1.asDiagonal()) * Omega.transpose() ;
  Eigen::MatrixXd M1rst = exp_mu_f0.asDiagonal() ;
  Eigen::MatrixXd Mbetw = (exp_mu_f1_f0.asDiagonal()) * Omega.transpose() ;
  
  // Compute macro stationary distribution:
  Eigen::MatrixXd Omega_2h = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd stat_distri = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Omega_2h = Omega.transpose() ;
  for (int i = 0; i < 10; i++){
    Omega_2h = Omega_2h * Omega_2h ;
  }
  stat_distri = Omega_2h.col(0) ;
  
  Eigen::MatrixXd all_LT_Bth          = Eigen::MatrixXd::Zero(nb_m,maxH) ;
  Eigen::MatrixXd all_LT_rth          = Eigen::MatrixXd::Zero(nb_m,maxH) ;
  Eigen::MatrixXd all_LT_ExpReturn_th = Eigen::MatrixXd::Zero(nb_m,maxH) ;
  
  Eigen::MatrixXd Bth          = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd rth          = Eigen::MatrixXd::Zero(nb_m,1) ;
  Eigen::MatrixXd ExpReturn_th = Eigen::MatrixXd::Zero(nb_m,1) ;
  
  // To compute expected returns:
  Eigen::MatrixXd exp_index   = Eigen::MatrixXd::Zero(nb_m,1) ;
  exp_index                   = ((mult(mu_pi,kappa_pi)+mult(mu_y,kappa_y)).array()).exp() ;
  Eigen::MatrixXd D_exp_indextOmega   = (exp_index.asDiagonal()) * Omega.transpose() ;
  Eigen::MatrixXd D_exp_indextOmega_h = Eigen::MatrixXd::Identity(nb_m,nb_m) ;
  
  // Compute bond prices:
  Eigen::MatrixXd Mbetw_h = Eigen::MatrixXd::Identity(nb_m, nb_m) ;
  for (int h = 0; h < maxH; h++){
    Bth = vec_1_m.transpose() * Mlast * Mbetw_h * M1rst ;
    Mbetw_h = Mbetw_h * Mbetw ; // for next iteration
    rth = - (Bth.array()).log() ;
    rth = mult(rth,1.0/(h+1)) ;
    
    D_exp_indextOmega_h = D_exp_indextOmega_h * D_exp_indextOmega ;
    ExpReturn_th = ((vec_1_m.transpose() * D_exp_indextOmega_h).transpose()).array() *
      (Bth.cwiseInverse()).array() ;
    
    all_LT_Bth.col(h) = Bth ;
    all_LT_rth.col(h) = rth ;
    all_LT_ExpReturn_th.col(h) = ExpReturn_th ;
  }
  
  return List::create(Named("all_LT_Bth")          = all_LT_Bth,
                      Named("all_LT_rth")          = all_LT_rth,
                      Named("all_LT_ExpReturn_th") = all_LT_ExpReturn_th,
                      Named("stat_distri")         = stat_distri) ;
}


// [[Rcpp::export]]
Rcpp::List compute_bond_prices(const Rcpp::List Model,
                               const int maxH,
                               const Eigen::MatrixXd indicators_x,
                               const Eigen::MatrixXd all_proba_def,
                               const Eigen::MatrixXd Probas){
  
  Rcpp::List res_LTprices    = compute_LTRF_bond_prices(Model, maxH) ;
  Eigen::MatrixXd all_LT_Bth = res_LTprices("all_LT_Bth") ;
  
  int nb_m         = all_LT_Bth.rows() ; // number of macro regimes
  int nb_states    = indicators_x.rows() ; // total number of states
  int nb_eps_m     = indicators_x.cols() ; // includes also macroeconomic regimes
  int nb_not_macro = 1.0 * nb_states / nb_m ; // number of regimes -- excluding macro regimes
  int nb_eps       = 1.0 * nb_eps_m / nb_m ;
  
  // compute nu associated with appropriate kappa's:
  double RR       = Model("RR") ;
  double gamma    = Model("gamma") ;
  double kappa_pi = Model("kappa_pi") ;
  double kappa_y  = Model("kappa_y") ;
  double nu_pi    = Model("nu_pi") ;
  double nu_y     = Model("nu_y") ;
  double nu = - gamma * nu_y + (kappa_pi - 1) * nu_pi + kappa_y * nu_y ;
  
  Eigen::MatrixXd vec_1_x     = Eigen::MatrixXd::Constant(nb_states,1,1) ;
  Eigen::MatrixXd vec_1_eps   = Eigen::MatrixXd::Constant(nb_eps,1,1) ;
  Eigen::MatrixXd vec_1_m_eps = Eigen::MatrixXd::Constant(nb_eps * nb_m,1,1) ;
  Eigen::MatrixXd vec_1_ddrr  = Eigen::MatrixXd::Constant(nb_not_macro, 1, 1) ;
  
  // extended version of all_LT_Bth (for all regimes, not only macro):
  Eigen::MatrixXd all_LT_Bth_ext = kronecker_cpp(all_LT_Bth,vec_1_ddrr) ; 
  
  Eigen::MatrixXd all_Bth       = Eigen::MatrixXd::Zero(nb_states,maxH) ;
  Eigen::MatrixXd all_rth       = Eigen::MatrixXd::Zero(nb_states,maxH) ;
  
  Eigen::MatrixXd all_prob_def  = Eigen::MatrixXd::Zero(nb_states,maxH) ;
  
  Eigen::MatrixXd p_h           = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd p_h_1         = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd all_p_h_1_tp1 = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  Eigen::MatrixXd p_h_aux       = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  
  Eigen::MatrixXd Bth              = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd rth              = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd Bth_aux          = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  Eigen::MatrixXd all_Bth_1_tp1    = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  Eigen::MatrixXd all_Bbarth_1_tp1 = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  
  Eigen::MatrixXd Bth_1            = Eigen::MatrixXd::Ones(nb_states,1) ;
  Eigen::MatrixXd B_bar_th_1       = Eigen::MatrixXd::Ones(nb_states,1) ;
  
  // s.d.f. specification (associated with composite index)
  Rcpp::List res_SDF = compute_SDF(Model) ;
  Eigen::MatrixXd mu_f0 = res_SDF("mu_f0") ;
  Eigen::MatrixXd mu_f1 = res_SDF("mu_f1") ;
  
  Eigen::MatrixXd all_f_tp1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  all_f_tp1 = vec_1_x * kronecker_cpp(mu_f1.transpose(), vec_1_eps.transpose()) +
    kronecker_cpp(mu_f0, vec_1_ddrr) * vec_1_m_eps.transpose() ;
  
  Eigen::MatrixXd all_exp_f_tp1     = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  Eigen::MatrixXd all_exp_f_tp1_nu  = Eigen::MatrixXd::Zero(nb_states,nb_eps * nb_m) ;
  all_exp_f_tp1    = (all_f_tp1.array()).exp() ;
  
  double RR_nu = RR * exp(nu) ;
  
  for (int h = 0; h < maxH; h++){
    all_p_h_1_tp1 = fill_from_indic(indicators_x, p_h_1) ;
    p_h_aux       =  all_p_h_1_tp1.array() + 
      add(- all_p_h_1_tp1,1).array() * all_proba_def.array() ;
    p_h_aux =  p_h_aux.array() * Probas.array() ;
    p_h = p_h_aux * vec_1_m_eps ;
    all_prob_def.col(h) = p_h ;
    p_h_1 = p_h ;
    
    all_Bth_1_tp1    = fill_from_indic(indicators_x, Bth_1) ;
    all_Bbarth_1_tp1 = fill_from_indic(indicators_x, B_bar_th_1) ;
    
    Bth_aux = (mult(all_Bbarth_1_tp1,RR_nu).array() - all_Bth_1_tp1.array()).array() * 
      all_proba_def.array() ;
    Bth_aux = all_exp_f_tp1.array() * Bth_aux.array() ;
    Bth_aux = all_exp_f_tp1.array() * all_Bth_1_tp1.array() + Bth_aux.array() ;
    Bth_aux =  Bth_aux.array() * Probas.array() ;
    Bth = Bth_aux * vec_1_m_eps ;
    
    rth = - (Bth.array()).log() ;
    rth = mult(rth,1.0/(h+1)) ;
    
    all_Bth.col(h) = Bth ;
    all_rth.col(h) = rth ;
    
    // for next iteration:
    Bth_1 = Bth ;
    B_bar_th_1 = all_LT_Bth_ext.col(h) ;
  }
  
  return List::create(Named("all_Bth")       = all_Bth,
                      Named("all_rth")       = all_rth,
                      Named("all_prob_def")  = all_prob_def,
                      Named("all_exp_f_tp1") = all_exp_f_tp1) ;
}




// [[Rcpp::export]]
Eigen::MatrixXd compute_stat_distri(const Rcpp::List Model){
  
  Eigen::MatrixXd Omega  = Model("Omega") ;
  int nb_m               = Omega.rows() ;
  
  Eigen::MatrixXd vec_1_m     = Eigen::MatrixXd::Constant(nb_m, 1, 1) ;
  
  // Compute macro stationary distribution:
  Eigen::MatrixXd Omega_2h = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd stat_distri = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Omega_2h = Omega.transpose() ;
  for (int i = 0; i < 10; i++){
    Omega_2h = Omega_2h * Omega_2h ;
  }
  stat_distri = Omega_2h.col(0) ;
  
  return stat_distri ;
}


// [[Rcpp::export]]
Rcpp::List compute_stat_distri_and_rbar(const Rcpp::List Model){
  
  Eigen::MatrixXd Omega  = Model("Omega") ;
  
  double chi       = Model("chi") ;
  
  int nb_m = Omega.rows() ;
  
  Eigen::MatrixXd vec_1_m     = Eigen::MatrixXd::Constant(nb_m, 1, 1) ;
  
  // Compute macro stationary distribution:
  Eigen::MatrixXd stat_distri = compute_stat_distri(Model) ;
  
  Eigen::MatrixXd expmuf0  = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd expmuf1  = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd expmuf01 = Eigen::MatrixXd::Zero(nb_m, 1) ;
  
  // s.d.f. specification (associated with composite index)
  Rcpp::List res_SDF = compute_SDF(Model) ;
  Eigen::MatrixXd mu_f0      = res_SDF("mu_f0") ;
  Eigen::MatrixXd mu_f1      = res_SDF("mu_f1") ;
  
  expmuf0  = (mu_f0.array()).exp() ;
  expmuf1  = (mu_f1.array()).exp() ;
  expmuf01 = ((mu_f0 + mu_f1).array()).exp() ;
  
  Eigen::MatrixXd Dlast = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Dbetw = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd D1rst = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Mlast = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Mbetw = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd M1rst = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  
  Dlast = expmuf1.asDiagonal() ;
  Dbetw = expmuf01.asDiagonal() ;
  D1rst = expmuf0.asDiagonal() ;
  Mlast = Dlast * Omega.transpose() ;
  Mbetw = Dbetw * Omega.transpose() ;
  M1rst = D1rst ;
  
  Eigen::MatrixXd I_m = Eigen::MatrixXd::Identity(nb_m, nb_m) ;
  Eigen::MatrixXd OnepChiPstar = Eigen::MatrixXd::Zero(1, nb_m) ;
  Eigen::MatrixXd Pstar = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd rstar = Eigen::MatrixXd::Zero(nb_m, 1) ;
  
  Pstar = (vec_1_m.transpose() * Mlast * (I_m - mult(Mbetw,chi)).inverse()  * M1rst).transpose() ;
  rstar = add(Pstar.cwiseInverse(), - 1 + chi) ;
  OnepChiPstar = add(mult(Pstar,chi),1) ;
  
  double r_bar = 1/(Pstar * stat_distri)(0,0) - 1 + chi ;
  

  return List::create(Named("r_bar")        = r_bar,
                      Named("stat_distri")  = stat_distri,
                      Named("OnepChiPstar") = OnepChiPstar,
                      Named("Pstar")        = Pstar,
                      Named("rstar")        = rstar
  ) ;
}


// [[Rcpp::export]]
Rcpp::List solve_ToyModel(const Eigen::MatrixXd all_d,
                          const Eigen::MatrixXd all_rr,
                          const Eigen::MatrixXd all_eps,
                          const Eigen::MatrixXd proba_eps,
                          const Rcpp::List Model,
                          const int nb_iter
){
  
  Eigen::MatrixXd Omega  = Model("Omega") ;
  Eigen::MatrixXd mu_pi  = Model("mu_pi") ;
  Eigen::MatrixXd mu_y   = Model("mu_y") ;
  Eigen::MatrixXd mu_eta = Model("mu_eta") ;
  
  double gamma    = Model("gamma") ;
  double chi      = Model("chi") ;
  double beta     = Model("beta") ;
  double d_star   = Model("d_star") ;
  double alpha    = Model("alpha") ;
  double s_star   = Model("s_star") ;
  double RR       = Model("RR") ;
  double kappa_pi = Model("kappa_pi") ;
  double kappa_y  = Model("kappa_y") ;
  double nu_pi    = Model("nu_pi") ;
  double nu_y     = Model("nu_y") ;
  
  int nb_grid_d  = all_d.rows() ;
  int nb_grid_rr = all_rr.rows() ;
  int nb_m       = Omega.rows() ;
  int nb_states  = nb_grid_d * nb_grid_d * nb_grid_rr * nb_m ;
  int nb_eps     = all_eps.rows() ;
  
  //XXXXX
  //Rcpp::List res_stat_distr_and_rbar = compute_stat_distri_and_rbar(Model) ;
  Eigen::MatrixXd expmuf0  = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd expmuf1  = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd expmuf01 = Eigen::MatrixXd::Zero(nb_m, 1) ;
  
  Eigen::MatrixXd Dlast = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Dbetw = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd D1rst = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Mlast = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd Mbetw = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  Eigen::MatrixXd M1rst = Eigen::MatrixXd::Zero(nb_m, nb_m) ;
  
  Eigen::MatrixXd vec_1_d     = Eigen::MatrixXd::Constant(nb_grid_d, 1, 1) ;
  Eigen::MatrixXd vec_1_rr    = Eigen::MatrixXd::Constant(nb_grid_rr, 1, 1) ;
  Eigen::MatrixXd vec_1_eps   = Eigen::MatrixXd::Constant(nb_eps, 1, 1) ;
  Eigen::MatrixXd vec_1_m     = Eigen::MatrixXd::Constant(nb_m, 1, 1) ;
  Eigen::MatrixXd vec_1_m_eps = Eigen::MatrixXd::Constant(nb_m*nb_eps, 1, 1) ;
  Eigen::MatrixXd vec_1_d2r   = Eigen::MatrixXd::Constant(nb_grid_d*nb_grid_d*nb_grid_rr, 1, 1) ;
  Eigen::MatrixXd vec_1_drm   = Eigen::MatrixXd::Constant(nb_grid_d*nb_grid_rr*nb_m, 1, 1) ;
  Eigen::MatrixXd vec_1_x     = Eigen::MatrixXd::Constant(nb_states, 1, 1) ;
  
  Eigen::MatrixXd d   = Eigen::MatrixXd::Zero(nb_states, 1) ;
  Eigen::MatrixXd d_1 = Eigen::MatrixXd::Zero(nb_states, 1) ;
  Eigen::MatrixXd rr  = Eigen::MatrixXd::Zero(nb_states, 1) ;
  Eigen::MatrixXd Pi  = Eigen::MatrixXd::Zero(nb_states, 1) ;
  Eigen::MatrixXd Dy  = Eigen::MatrixXd::Zero(nb_states, 1) ;
  Eigen::MatrixXd s_m = Eigen::MatrixXd::Zero(nb_states, 1) ;
  
  Eigen::MatrixXd all_rr_tp1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_d_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  
  Eigen::MatrixXd indicators_d_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd indicators_m_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd indicators_d_t    = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd indicators_rr_tp1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd indicators_x      = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_lambdas       = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_proba_def     = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_q_tp1         = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_q0_tp1        = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd E                 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd Q                 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  
  Eigen::MatrixXd stat_distri = compute_stat_distri(Model) ;
  
  // s.d.f. specification (associated with composite index)
  Rcpp::List res_SDF = compute_SDF(Model) ;
  Eigen::MatrixXd mu_f0 = res_SDF("mu_f0") ;
  Eigen::MatrixXd mu_f1 = res_SDF("mu_f1") ;
  
  expmuf0  = (mu_f0.array()).exp() ;
  expmuf1  = (mu_f1.array()).exp() ;
  expmuf01 = ((mu_f0 + mu_f1).array()).exp() ;
  
  Dlast = expmuf1.asDiagonal() ;
  Dbetw = expmuf01.asDiagonal() ;
  D1rst = expmuf0.asDiagonal() ;
  Mlast = Dlast * Omega.transpose() ;
  Mbetw = Dbetw * Omega.transpose() ;
  M1rst = D1rst ;
  
  Eigen::MatrixXd I_m = Eigen::MatrixXd::Identity(nb_m, nb_m) ;
  Eigen::MatrixXd Pstar        = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd rstar        = Eigen::MatrixXd::Zero(nb_m, 1) ;
  Eigen::MatrixXd OnepChiPstar = Eigen::MatrixXd::Zero(1, nb_m) ;
  Pstar        = (vec_1_m.transpose() * Mlast * (I_m - mult(Mbetw,chi)).inverse() * M1rst).transpose() ;
  OnepChiPstar = add(mult(Pstar.transpose(),chi),1) ;
  //XXXXX
  rstar = add(Pstar.cwiseInverse(), - 1 + chi) ; // will be used to initialize q

  Eigen::MatrixXd Mat_1 = Eigen::MatrixXd::Constant(nb_states, nb_m * nb_eps, 1) ;
  
  d   = kronecker_cpp(kronecker_cpp(kronecker_cpp(vec_1_m,vec_1_rr),vec_1_d),all_d) ;
  d_1 = kronecker_cpp(kronecker_cpp(kronecker_cpp(vec_1_m,vec_1_rr),all_d)  ,vec_1_d) ;
  rr  = kronecker_cpp(kronecker_cpp(kronecker_cpp(vec_1_m,all_rr)  ,vec_1_d),vec_1_d) ;
  Pi  = kronecker_cpp(kronecker_cpp(kronecker_cpp(mu_pi,  vec_1_rr),vec_1_d),vec_1_d) ;
  Dy  = kronecker_cpp(kronecker_cpp(kronecker_cpp(mu_y,   vec_1_rr),vec_1_d),vec_1_d) ;
  s_m = kronecker_cpp(kronecker_cpp(kronecker_cpp(mu_eta, vec_1_rr),vec_1_d),vec_1_d) ;
  
  double nu = - gamma * nu_y + (kappa_pi - 1) * nu_pi + kappa_y * nu_y ;
  
  Eigen::MatrixXd all_eps_tp1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_Pi_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_Dy_tp1  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_f_tp1   = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_OnepChiPstar = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  
  // Eigen::MatrixXd vec_1_m_transp = vec_1_m.transpose() ;
  // Eigen::MatrixXd all_eps_transp = all_eps.transpose() ;
  
  all_eps_tp1 = vec_1_x * kronecker_cpp(vec_1_m.transpose(), all_eps.transpose()) ;
  all_Pi_tp1  = vec_1_x * kronecker_cpp(mu_pi.transpose(),   vec_1_eps.transpose()) ;
  all_Dy_tp1  = vec_1_x * kronecker_cpp(mu_y.transpose(),    vec_1_eps.transpose()) ;
  all_f_tp1   = vec_1_x * kronecker_cpp(mu_f1.transpose(),   vec_1_eps.transpose()) +
    kronecker_cpp(kronecker_cpp(kronecker_cpp(mu_f0, vec_1_rr), vec_1_d),vec_1_d) * vec_1_m_eps.transpose() ;
  all_OnepChiPstar = vec_1_x * kronecker_cpp(OnepChiPstar, vec_1_eps.transpose()) ;
  
  Eigen::MatrixXd all_d_t   = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_d_t_1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_rr_t  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_Pi_t  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_Dy_t  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_q_t   = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_q0_t  = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Eigen::MatrixXd all_s_m_t = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  
  Eigen::MatrixXd Matrix_beta = Eigen::MatrixXd::Constant(nb_states, nb_m * nb_eps, beta) ;
  
  Eigen::MatrixXd q     = Eigen::MatrixXd::Zero(nb_states, 1) ; // vector of sovereign yields
  Eigen::MatrixXd q0    = Eigen::MatrixXd::Zero(nb_states, 1) ; // risk-free yields
  
  q  = kronecker_cpp(kronecker_cpp(kronecker_cpp(add(rstar,.001), vec_1_rr),vec_1_d),vec_1_d) ;
  q0 = kronecker_cpp(kronecker_cpp(kronecker_cpp(add(rstar,.000), vec_1_rr),vec_1_d),vec_1_d) ;
  
  Eigen::MatrixXd AUX = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  
  all_d_t   = d   * vec_1_m_eps.transpose() ;
  all_d_t_1 = d_1 * vec_1_m_eps.transpose() ;
  all_rr_t  = rr  * vec_1_m_eps.transpose() ;
  all_Pi_t  = Pi  * vec_1_m_eps.transpose() ;
  all_Dy_t  = Dy  * vec_1_m_eps.transpose() ;
  all_q_t   = q   * vec_1_m_eps.transpose() ;
  all_q0_t  = q0  * vec_1_m_eps.transpose() ;
  all_s_m_t = s_m * vec_1_m_eps.transpose() ;
  
  Eigen::MatrixXd Probas   = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  Probas = kronecker_cpp(kronecker_cpp(Omega,proba_eps.transpose()),vec_1_d2r) ;
  
  Eigen::MatrixXd q_iter_1 = q ; // used to check convergence over iterations
  Eigen::MatrixXd q_chge   = Eigen::MatrixXd::Zero(nb_states, 1) ; // used to check convergence
  
  Eigen::MatrixXd all_zeta_t   = ((mult(all_Pi_t,  kappa_pi - 1) + mult(all_Dy_t,  kappa_y - 1)).array()).exp() ;
  Eigen::MatrixXd all_zeta_tp1 = ((mult(all_Pi_tp1,kappa_pi - 1) + mult(all_Dy_tp1,kappa_y - 1)).array()).exp() ;
  
  Eigen::MatrixXd all_eta_tp1 = Eigen::MatrixXd::Zero(nb_states, nb_m * nb_eps) ;
  all_eta_tp1 = all_eps_tp1 +
    vec_1_x * kronecker_cpp(mu_eta.transpose(), vec_1_eps.transpose()) - all_s_m_t ;
  
  Eigen::MatrixXd seq1d = seqEigen(1,nb_grid_d) ;
  Eigen::MatrixXd seq1m = seqEigen(1,nb_m) ;
  
  indicators_d_t   = kronecker_cpp(vec_1_drm,seq1d) * vec_1_m_eps.transpose() ;
  indicators_m_tp1 = vec_1_x * kronecker_cpp(seq1m.transpose(),vec_1_eps.transpose()) ;
  
  if(nb_iter>0){
    
    for (int i = 0; i < nb_iter; i++){
      
      all_q_t = q * vec_1_m_eps.transpose() ;
      
      all_rr_tp1 = all_q_t.array() * all_zeta_tp1.array() *
        (all_d_t - mult(all_zeta_t.array() * all_d_t_1.array(),chi)).array() +
        mult(all_zeta_tp1.array() * all_rr_t.array(),chi).array() ;
      
      all_d_tp1  = all_zeta_tp1.array() * all_d_t.array() -
        mult(add(all_d_t, - d_star),beta).array() -
        all_eta_tp1.array() + all_rr_tp1.array() ;
      
      // all_d_tp1  = all_zeta_tp1.array() * all_d_t.array() -
      //   Matrix_beta.array() - all_eta_tp1.array() + all_rr_tp1.array() ;
      
      indicators_d_tp1  = find_closest_in_vec(all_d_tp1,all_d) ;
      indicators_rr_tp1 = find_closest_in_vec(all_rr_tp1,all_rr) ;
      
      indicators_x = indicators_d_tp1.array() +
        mult(indicators_d_t - Mat_1,nb_grid_d).array() +
        mult(indicators_rr_tp1 - Mat_1,nb_grid_d*nb_grid_d).array() +
        mult(indicators_m_tp1 - Mat_1,nb_grid_d*nb_grid_d*nb_grid_rr).array() ;
      
      // Compute probabilities of default:
      all_lambdas = pmax_cpp(
        add(mult(add(all_d_t, - d_star),beta).array() + all_eta_tp1.array(), - s_star),0) ;
      // all_lambdas = pmax_cpp(add(all_d_t, - d_star),0) ;
      
      all_proba_def = add(- ((mult(all_lambdas,-alpha)).array()).exp(),1) ;
      
      // Update q (sovereign yields):
      all_q_tp1 = fill_from_indic(indicators_x, q) ;
      
      AUX = add(all_q_tp1,1).array() * (add(all_q_tp1,1 - chi).cwiseInverse()).array() ;
      
      E = ((all_f_tp1.array()).exp()).array() * (
        add(all_q_tp1,1).array() * (add(all_q_tp1,1 - chi).cwiseInverse()).array() +
          all_proba_def.array() * (mult(all_OnepChiPstar,exp(nu)*RR).array() -
          add(all_q_tp1,1).array() * (add(all_q_tp1,1 - chi).cwiseInverse()).array()).array()
      ).array() ;
      
      Q = E.array() * Probas.array() ;
      q = add((Q * vec_1_m_eps).cwiseInverse(),chi - 1) ; // update q

      // Update q0 (risk-free yield, obtained for RR=1):
      all_q0_tp1 = fill_from_indic(indicators_x, q0) ;
      
      E = ((all_f_tp1.array()).exp()).array() * (
        add(all_q0_tp1,1).array() * (add(all_q0_tp1,1 - chi).cwiseInverse()).array() +
          all_proba_def.array() * (mult(all_OnepChiPstar,exp(nu)).array() -
          add(all_q0_tp1,1).array() * (add(all_q0_tp1,1 - chi).cwiseInverse()).array()).array()
      ).array() ;
      
      Q = add(E.cwiseInverse(),chi - 1).array() * Probas.array() ;
      q0 = Q * vec_1_m_eps ; // update q0
      
      q_chge = q - q_iter_1 ;
      q_iter_1 = q; // for next iteration
    }
  }
  return List::create(Named("q")  = q,
                      Named("q0") = q0,
                      Named("d") = d,
                      Named("d_1") = d_1,
                      Named("rr") = rr,
                      Named("Pi") = Pi,
                      Named("Dy") = Dy,
                      Named("q_chge") = q_chge,
                      Named("indicators_x") = indicators_x,
                      Named("all_f_tp1") = all_f_tp1,
                      Named("mu_f0") = mu_f0,
                      Named("mu_f1") = mu_f1,
                      Named("nu")   = nu,
                      Named("all_proba_def") = all_proba_def,
                      Named("Probas") = Probas,
                      Named("stat_distri") = stat_distri,
                      Named("OnepChiPstar") = OnepChiPstar,
                      Named("Pstar") = Pstar,
                      Named("rstar") = rstar
  ) ;
}


// [[Rcpp::export]]
Eigen::MatrixXd compute_proba_def(const int maxH,
                                  const Eigen::MatrixXd indicators_x,
                                  const Eigen::MatrixXd all_proba_def,
                                  const Eigen::MatrixXd Probas){
  
  int nb_states = indicators_x.rows() ;
  int nb_eps    = indicators_x.cols() ;
  
  Eigen::MatrixXd vec_1_eps = Eigen::MatrixXd::Constant(nb_eps,1,1) ;
  
  Eigen::MatrixXd all_prob_def  = Eigen::MatrixXd::Zero(nb_states,maxH) ;
  Eigen::MatrixXd p_h           = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd p_h_1         = Eigen::MatrixXd::Zero(nb_states,1) ;
  Eigen::MatrixXd all_p_h_1_tp1 = Eigen::MatrixXd::Zero(nb_states,nb_eps) ;
  Eigen::MatrixXd p_h_aux       = Eigen::MatrixXd::Zero(nb_states,nb_eps) ;
  
  for (int h = 0; h < maxH; h++){
    all_p_h_1_tp1    = fill_from_indic(indicators_x, p_h_1) ;
    p_h_aux          =  all_p_h_1_tp1.array() + 
      add(- all_p_h_1_tp1,1).array() * all_proba_def.array() ;
    p_h_aux =  p_h_aux.array() * Probas.array() ;
    p_h = p_h_aux * vec_1_eps ;
    all_prob_def.col(h) = p_h ;
    p_h_1 = p_h ;
  }
  
  return all_prob_def ;
}




// [[Rcpp::export]]
Eigen::MatrixXd compute_multivariate_normal(const Eigen::MatrixXd epsilon_matrix,
                                            const Eigen::MatrixXd Covariance){
  
  int nb_dates = epsilon_matrix.rows() ;
  int nb_eps   = epsilon_matrix.cols() ;
  int nb_eps2  = nb_eps * nb_eps ;
  
  Eigen::MatrixXd vec_1_dates = Eigen::MatrixXd::Constant(nb_dates,1,1) ;
  Eigen::MatrixXd vec_1_eps   = Eigen::MatrixXd::Constant(nb_eps,1,1) ;
  Eigen::MatrixXd vec_1_eps2  = Eigen::MatrixXd::Constant(nb_eps2,1,1) ;
  
  Eigen::MatrixXd loglik = Eigen::MatrixXd::Zero(nb_dates,1) ;
  
  double det_Covariance = Covariance.determinant() ;
  
  Eigen::MatrixXd Covariance_inv = Covariance.inverse() ;
  
  Eigen::MatrixXd vec_Covariance_inv = Eigen::MatrixXd::Zero(nb_dates,nb_eps2) ;
  vec_Covariance_inv = vec_1_dates * Vec(Covariance_inv).transpose() ;
  
  vec_Covariance_inv = kronecker_cpp(epsilon_matrix,vec_1_eps.transpose()).array() *
    vec_Covariance_inv.array() * kronecker_cpp(vec_1_eps.transpose(),epsilon_matrix).array() ;
  
  double pi = 3.141592653589793238462643383 ;
  
  double aux = -.5 * nb_eps * log(2*pi) -.5 * log(det_Covariance) ;
  
  loglik = add(mult(vec_Covariance_inv * vec_1_eps2, -0.5),aux) ;
  
  return loglik ;
}

// [[Rcpp::export]]
Rcpp::List KH_filter(const Eigen::MatrixXd F,
                     const Eigen::MatrixXd M,
                     const Eigen::MatrixXd N,
                     const Eigen::MatrixXd Omega){
  
  // Kitagawa-Hamilton filter --------------------------------------------------
  // The model is:
  // F_t = M·z_t + diag(N·z_t)·epsilon_t
  // where:
  //   * z_t is a selection vector following a homogenous Markov_switching process
  //     with a matrix of transition probabilities Omega (lines sum to one), and
  //   * epsilon_t is i.i.d. N(0,Id)
  // ---------------------------------------------------------------------------
  
  int J          = M.cols() ; // number of regimes
  int nb_dates   = F.rows() ;
  int nb_obsvar  = F.cols() ;
  
  Eigen::MatrixXd ksi_matrix   = Eigen::MatrixXd::Zero(nb_dates,J) ;
  Eigen::MatrixXd ksi_1_matrix = Eigen::MatrixXd::Zero(nb_dates,J) ;
  Eigen::MatrixXd eta_matrix   = Eigen::MatrixXd::Zero(nb_dates,J) ;
  Eigen::MatrixXd ksi_t        = Eigen::MatrixXd::Zero(J,1) ;
  Eigen::MatrixXd eta_t        = Eigen::MatrixXd::Zero(J,1) ;
  
  Eigen::MatrixXd vec_1_dates = Eigen::MatrixXd::Constant(nb_dates,1,1) ;
  Eigen::MatrixXd vec_1_J     = Eigen::MatrixXd::Constant(J,1,1) ;
  
  // Compute log-lik conditional on each regime:
  Eigen::MatrixXd variances4regime = Eigen::MatrixXd::Zero(nb_obsvar,1) ;
  Eigen::MatrixXd Covariance       = Eigen::MatrixXd::Zero(nb_obsvar,nb_obsvar) ;
  for (int j = 0; j < J; j++){
    variances4regime = (N.col(j)).array() * (N.col(j)).array() ;
    Covariance = variances4regime.asDiagonal() ;
    eta_matrix.col(j) = compute_multivariate_normal(F - vec_1_dates * (M.col(j)).transpose(),Covariance) ;
  }
  
  // Compute stationary distribution:
  Eigen::MatrixXd Omega_2h = Eigen::MatrixXd::Zero(J,J) ;
  Eigen::MatrixXd stat_distri = Eigen::MatrixXd::Zero(J, 1) ;
  Omega_2h = Omega.transpose() ;
  for (int i = 0; i < 10; i++){
    Omega_2h = Omega_2h * Omega_2h ;
  }
  stat_distri = Omega_2h.col(0) ;
  ksi_t = stat_distri ;
  
  for (int t = 0; t < nb_dates; t++){
    ksi_1_matrix.row(t) = ksi_t.transpose() ; // will be used to compute log-likelihood
    
    eta_t = (eta_matrix.row(t)).transpose() ; // eta_matrix is LOG-likelihood
    eta_t = eta_t.array().exp() ;
    // Use update formula:
    ksi_t = (Omega.transpose() * ksi_t).array() * eta_t.array() ;
    double normalisation_factor = (vec_1_J.transpose() * ksi_t)(0,0) ;
    ksi_t = mult(ksi_t,1/normalisation_factor) ;
    
    ksi_matrix.row(t) = ksi_t.transpose() ;
    eta_matrix.row(t) = eta_t.transpose() ;
  }
  
  // Compute log-likelihood :
  Eigen::MatrixXd loglik_mat = Eigen::MatrixXd::Zero(nb_dates,J) ;
  Eigen::MatrixXd loglik_vec = Eigen::MatrixXd::Zero(nb_dates,1) ;
  loglik_mat = ((ksi_1_matrix * Omega).array() * eta_matrix.array()) ;
  loglik_vec = loglik_mat * vec_1_J ;
  loglik_vec = (loglik_vec.array()).log() ;
  double loglik = (loglik_vec.transpose() * vec_1_dates)(0,0) ;
  
  return List::create(Named("ksi_matrix") = ksi_matrix,
                      Named("eta_matrix") = eta_matrix,
                      Named("loglik_vec") = loglik_vec,
                      Named("loglik")     = loglik) ;
}



// [[Rcpp::export]]
Eigen::MatrixXd compute_uncond_distri(const Eigen::MatrixXd & indicators_x,
                                      const Eigen::MatrixXd & Probas,
                                      const int nb_iter){
  int nb_states = indicators_x.rows() ;
  int nb_eps_m  = indicators_x.cols() ;
  
  double p_ini = 1/(1.0 * nb_states) ;
  
  Eigen::MatrixXd p_1 = Eigen::MatrixXd::Constant(nb_states,1,p_ini) ;
  Eigen::MatrixXd p   = Eigen::MatrixXd::Zero(nb_states,1) ;
  
  for(int iter = 0; iter < nb_iter; iter++){
    p = mult(p,0) ;
    for(int i = 0; i < nb_states; i++){
      if(p_1(i,0) > 0){
        for(int j = 0; j < nb_eps_m; j++){
          int k = indicators_x(i,j) - 1;
          p(k,0) = p(k,0) + Probas(i,j) * p_1(i,0) ;
        }
      }
    }
    p_1 = p ;
  }
  return p ;
}

// [[Rcpp::export]]
Eigen::MatrixXd compute_distri_x(const Eigen::MatrixXd & all_x,
                                 const Eigen::MatrixXd & x,
                                 const Eigen::MatrixXd & p){
  int nb_states = x.rows() ;
  int nb_x      = all_x.rows() ;
  
  Eigen::MatrixXd px           = Eigen::MatrixXd::Zero(nb_x,1) ;
  Eigen::MatrixXd indicators_x = Eigen::MatrixXd::Zero(nb_states,1) ;
  
  indicators_x = find_closest_in_vec(x,all_x) ;
  
  for(int i = 0; i < nb_states; i++){
    px(indicators_x(i,0)-1,0) = px(indicators_x(i,0)-1,0) + p(i,0) ;
  }
  
  return px ;
}