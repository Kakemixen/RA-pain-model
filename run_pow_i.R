source("./general.R")
library(ggplot2)
library(rstan)

model_name = "pow_i"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 4

# get pars vector
paramList = c("RiskAversion","PainAvoidance","PainRetention","tau","log_lik","PredictedResponse")
dataList = get_dataList()

#output = sample_model(model_name, dataList, paramList, iterations, warmups, chains, init= list(
#                          chain_1 = list(
#                               "PainRetention"       = rep(1, 35)
#                                ),
#                          chain_2 = list(
#                               "PainRetention"       = rep(3.1, 35)
#                                ),
#                          chain_3 = list(
#                               "PainRetention"       = rep(2, 35)
#                                ),
#                          chain_4 = list(
#                               "PainRetention"       = rep(0.5, 35)
#                                )
#                      )
#)

BIC(output, dataList, 6)
PPC(output, dataList)
LOOIC(output)

## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("RiskAversion"))
traceplot(output, pars=c("PainAvoidance"))
traceplot(output, pars=c("PainRetention"))
traceplot(output, pars=c("tau"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("RiskAversion"))
stan_plot(output, pars=c("PainAvoidance"))
stan_plot(output, pars=c("PainRetention"))
stan_plot(output, pars=c("tau"))

