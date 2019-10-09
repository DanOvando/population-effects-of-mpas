
fit_depletion <- function(target_depletion,scientific_name, target_catch, linf,sim_years,burn_years,
                   num_patches,
                   larval_movement,
                   adult_movement,
                   density_dependence_form,
                   steepness)
{
  nlminb(
    c(log(1000*target_catch), log(200)),
    tune_f,
    target_depletion = target_depletion,
    target_catch = target_catch,
    scientific_name = scientific_name,
    linf = NA,
    alpha = 0.5,
    sim_years = sim_years,
    burn_years = burn_years,
    num_patches = num_patches,
    adult_movement = adult_movement,
    larval_movement = larval_movement,
    steepness = steepness,
    density_dependence_form = density_dependence_form
  )$par

}
