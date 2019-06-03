source("./general.R")
library(ggplot2)
library(rstan)

model_name = "const_h"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 1

# get pars vector
paramList = c("RiskAversion","PainAvoidance","tau","mu_p", "sigma_p", "mu_RiskAversion","mu_PainAvoidance","mu_tau","log_lik","PredictedResponse")
dataList = get_dataList()

output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

PPC(output, dataList, iterations-warmups)

## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("RiskAversion"))
traceplot(output, pars=c("PainAvoidance"))
traceplot(output, pars=c("tau"))
traceplot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("RiskAversion"))
stan_plot(output, pars=c("PainAvoidance"))
# plot(density(extracted$PainAvoidance[,,3]), col="red")
stan_plot(output, pars=c("tau"))
stan_plot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))
stan_dens(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))
# stan_plot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))
# stan_dens(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))

