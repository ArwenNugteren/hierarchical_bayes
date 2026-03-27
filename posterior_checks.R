stanFit <- stanfit_timestep_216

cat("Max Rhat:  ", max(summary(stanFit)$summary[, "Rhat"],  na.rm = TRUE), "\n")
cat("Min n_eff: ", min(summary(stanFit)$summary[, "n_eff"], na.rm = TRUE), "\n")

View(summary(stanFit)$summary)


# hyperparameters
pairs(stanFit, pars = c("mu_alpha", "mu_gamma"))

# lower level parameters
# get all parameter names from the fit
all_pars <- names(stanFit)

# grab all alpha-related parameters
alpha_pars <- all_pars[grepl("alpha_state", all_pars)]
pairs(stanFit, pars = alpha_pars)

alpha_means <- summary(stanFit)$summary[alpha_pars, "mean"]
print(alpha_means)

# grab all gamma-related parameters
gamma_pars <- all_pars[grepl("gamma_state", all_pars)]
pairs(stanFit, pars = gamma_pars)

gamma_means <- summary(stanFit)$summary[gamma_pars, "mean"]
print(gamma_means)

# grab all l0-related parameters
l0_pars <- all_pars[grepl("l0", all_pars)]

pairs(stanFit, pars = l0_pars)