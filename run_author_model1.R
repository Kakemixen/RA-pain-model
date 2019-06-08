source("./general.R")
library(ggplot2)
library(rstan)

model_name = "author1"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 2

# get pars vector
paramList = c("beta_mu","theta", "ddb", "PointPosteriors", "PredictedResponse", "log_lik")
dataList = get_chris_dataList()


output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

#parameters <- rstan::extract(output)


chris_BIC(output, dataList, 1, iterations-warmups)
LOOIC(output)

## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("beta_mu"))
traceplot(output, pars=c("theta"))
traceplot(output, pars=c("ddb"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("beta_mu"))
stan_plot(output, pars=c("theta"))
stan_plot(output, pars=c("ddb"))

