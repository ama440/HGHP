//' Obtain horseshoe SD For levels of factors
//' @param deltalvl_raw raw form of delta
//' @param delta_grp    group shrinkage factor
//' @param lvlgrp_indx  integer array giving to which level each level belongs
//' @param nlevelparms  number of total level parameters
//' @return shrinkage   parameter for each level of each factor
vector get_delta_lvl(vector deltalvl_raw, vector delta_grp, array[] int lvlgrp_indx, int nlevelparms) {
  vector[nlevelparms] delta_lvl;
  for ( l in 1:nlevelparms )
    delta_lvl[l] = delta_grp[lvlgrp_indx[l]] * tan(deltalvl_raw[l]);
  return(delta_lvl);
}
//' Get standard deviation of all regression coefficients
//' @param delta_lvl     horseshoe SD for levels of factors
//' @param tau           global shrinkage parameters
//' @param lambda        local shrinkage parameters
//' @param lvl_indx      pxG integer array giving to which level of which group(s) each regression parameter belongs
//' @param indx_select   integer array giving which variables to perform variable selection
//' @param indx_force    integer array giving whdich variables are forced into the model
//' @param beta_force_sd standard deviation for normal prior on forced-in regression coefficients
//' @param p             number of regression coefficients
//' @param sd_force      vector giving standard deviation of forced-in regression coefficients
//' @param G             total number of groups
//' @return p-dim vector giving standard deviations for regression coefficients
vector get_beta_sd(
  vector delta_lvl, real tau, vector lambda
  , array[,] int lvl_indx, array[] int indx_select, array[] int indx_force
  , real beta_force_sd, int p, vector sd_force, int G
) {
  vector[p] sd_beta;
  // prior SD for forced-in betas
  sd_beta[indx_force] = sd_force;
  // for variable selection, initialize SD to tau * lambda_j
  sd_beta[indx_select] = tau * lambda;
  // multiply SD by level-level shrinkage factors
  for ( g in 1:G ) {
    for ( j in indx_select ) {
      if ( lvl_indx[j,g] > 0)
        sd_beta[j] *= delta_lvl[lvl_indx[j,g]];
    }
  }
  return(sd_beta);
}