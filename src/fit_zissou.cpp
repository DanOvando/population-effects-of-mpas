
// TMB attempt at fitting ahnold
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  /////////load data/////////

  // seen data

  DATA_MATRIX(x_seen_non_nested); // non nested part of data

  // DATA_MATRIX(x_seen_did); // difference in difference terms

  DATA_MATRIX(x_seen_year_species); // year random effects

  DATA_MATRIX(x_seen_region_cluster);

  DATA_IVECTOR(year_species_index); // index the same rows as x_seen showing what species is in that row

  DATA_IVECTOR(seen_species_index); // index the same rows as x_seen showing what species is in that row

  DATA_IVECTOR(region_cluster_index); // index the same rows as x_seen showing what species is in that row

  DATA_VECTOR(log_density); // observed log densities

  // seeing data

  DATA_MATRIX(x_seeing_non_nested); // non nested part of data

  DATA_MATRIX(x_seeing_year_species); // year random effects

  DATA_MATRIX(x_seeing_region_cluster); // REGION CLUSTERS

  DATA_VECTOR(any_seen); // observed log densities



  // standardized matrices data

  DATA_MATRIX(standard_non_nested); // standardized non nested part of data

  DATA_MATRIX(standard_year_species); // standardized did estimators without mpa

  DATA_MATRIX(standard_region_cluster); // STANDARD REGION CLUSTERS

  // did data

  DATA_MATRIX(non_nested_did_data);

  DATA_MATRIX(targeted_year_did_data);

  DATA_MATRIX(nontargeted_year_did_data);

  DATA_MATRIX(species_did_data);



  /////////define parameters/////////

  PARAMETER_VECTOR(seen_non_nested_betas); // NON DID BETAS

  // PARAMETER_VECTOR(seen_did_betas);  // DIFFERENCE IN DIFFERENCE BETAS

  PARAMETER_VECTOR(seen_density_species_sigma);

  PARAMETER_VECTOR(seen_year_species_betas);

  PARAMETER_VECTOR(seen_year_species_sigmas);

  PARAMETER_VECTOR(seen_region_cluster_betas);

  PARAMETER_VECTOR(seen_region_cluster_sigmas);


  // seeing parameters

  PARAMETER_VECTOR(seeing_non_nested_betas);

  PARAMETER_VECTOR(seeing_year_species_betas);

  PARAMETER_VECTOR(seeing_year_species_sigmas);

  PARAMETER_VECTOR(seeing_region_cluster_betas);

  PARAMETER_VECTOR(seeing_region_cluster_sigmas);

  // did parameters

  PARAMETER_VECTOR(non_nested_did_betas);

  PARAMETER_VECTOR(targeted_did_betas);

  PARAMETER_VECTOR(nontargeted_did_betas);

  PARAMETER_VECTOR(species_did_betas)

  PARAMETER(log_targeted_sigma);

  PARAMETER(log_nontargeted_sigma);

  PARAMETER(log_species_sigma);

  PARAMETER(log_did_sigma);

  /////////process parameters and data/////////

  Type nll; //blank storage for accumulated nll

  nll = 0;

  int i_max;
  /////////seen fish/////////

  matrix<Type> non_nested_effects = x_seen_non_nested * seen_non_nested_betas;

  // matrix<Type> did_effects = x_seen_did * seen_did_betas;

  matrix<Type> year_species_effects = x_seen_year_species * seen_year_species_betas;

  matrix<Type> region_cluster_effects = x_seen_region_cluster * seen_region_cluster_betas;

  matrix<Type> log_density_hat = non_nested_effects + year_species_effects + region_cluster_effects ; // + year_species_effects + region_cluster_effects;

  i_max = seen_year_species_betas.size();

  // std::cout << i_max << "\\n";

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(seen_year_species_betas(i), Type(0), exp(seen_year_species_sigmas(year_species_index(i) - 1)), true);


  } // close year species effects


  i_max = log_density.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(log_density(i), log_density_hat(i), exp(seen_density_species_sigma(seen_species_index(i) - 1)), true);

  } // close log density nll


  i_max = seen_region_cluster_betas.size();

  for (int i = 0; i < i_max ; i++){

    nll -= dnorm(seen_region_cluster_betas(i), Type(0), exp(seen_region_cluster_sigmas(region_cluster_index(i) - 1)), true);


  } // close region cluster effects


  /////////seeing fish/////////

  matrix<Type> seeing_non_nested_effects = x_seeing_non_nested * seeing_non_nested_betas;

  matrix<Type> seeing_year_species_effects = x_seeing_year_species * seeing_year_species_betas;

  matrix<Type> seeing_region_cluster_effects = x_seeing_region_cluster * seeing_region_cluster_betas;

  matrix<Type> logit_scale_prob_seeing = seeing_non_nested_effects + seeing_year_species_effects + seeing_region_cluster_effects; //+ seeing_year_species_effects + seeing_region_cluster_effects;

  vector<Type> prob_seeing =  1/ (1 + exp(-logit_scale_prob_seeing.array()));

  i_max = seeing_year_species_betas.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(seeing_year_species_betas(i), Type(0), exp(seeing_year_species_sigmas(year_species_index(i) - 1)), true);

  } // close year species effects

  nll -= sum(dbinom(any_seen,Type(1),prob_seeing, true));

  i_max = seeing_region_cluster_betas.size();

  for (int i = 0; i < i_max ; i++){

    nll -= dnorm(seeing_region_cluster_betas(i), Type(0), exp(seeing_region_cluster_sigmas(region_cluster_index(i) - 1)), true);


  } // close region cluster effects


  //// standardized abundances ////
  vector<Type> logit_standardized_yearly_prob_seeing = standard_non_nested * seeing_non_nested_betas + standard_year_species * seeing_year_species_betas  + standard_region_cluster * seeing_region_cluster_betas; //

  vector<Type> standardized_yearly_prob_seeing = 1 / (1 + exp(-logit_standardized_yearly_prob_seeing));

  vector<Type> log_standardized_abundance = standard_non_nested * seen_non_nested_betas + standard_year_species * seen_year_species_betas  + standard_region_cluster * seen_region_cluster_betas;

  vector<Type> standardized_abundance = exp(log_standardized_abundance);

  vector<Type> abundance_hat = standardized_yearly_prob_seeing * standardized_abundance;

  vector<Type> log_abundance_hat = log(abundance_hat);

  //// did model ////

  vector<Type> abundance_hat_hat = non_nested_did_data * non_nested_did_betas + targeted_year_did_data * targeted_did_betas + nontargeted_year_did_data * nontargeted_did_betas + species_did_data * species_did_betas;

  i_max = abundance_hat.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(log_abundance_hat(i), abundance_hat_hat(i), exp(log_did_sigma), true);

  } // close did thing


  i_max = targeted_did_betas.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(targeted_did_betas(i), Type(0), exp(log_targeted_sigma), true);

    nll -= dnorm(nontargeted_did_betas(i), Type(0), exp(log_nontargeted_sigma), true);

  } // close did thing


  i_max = species_did_betas.size();

  for (int i = 0; i < i_max; i++){

    nll -= dnorm(species_did_betas(i), Type(0), exp(log_species_sigma), true);

  } // close did thing


  vector<Type> mpa_effect = targeted_did_betas - nontargeted_did_betas;

  /////////outputs/////////


  REPORT(log_density_hat);

  REPORT(prob_seeing);

  REPORT(standardized_yearly_prob_seeing);

  REPORT(standardized_abundance);

  REPORT(abundance_hat);

  REPORT(mpa_effect);

  REPORT(abundance_hat_hat);

  ADREPORT(mpa_effect);

  ADREPORT(abundance_hat_hat);

  return nll;

}
