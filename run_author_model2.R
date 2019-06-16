source("./general.R")
library(ggplot2)
library(rstan)

model_name = "author2"
print("running model")
print(model_name)

iterations = 2000
warmups = 1000
chains = 4

# get pars vector
paramList = c("beta_mu", "beta", "beta_sig", "theta", "ddb", "PointPosteriors", "PredictedResponse", "log_lik")
dataList = get_chris_dataList()


output = sample_model(model_name, dataList, paramList, iterations, warmups, chains)

#parameters <- rstan::extract(output)
chris_BIC(output, dataList, 4)
LOOIC(output)



PPC_chris(output, dataList, iterations-warmups)

## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("beta_mu"))
traceplot(output, pars=c("beta"))
traceplot(output, pars=c("beta_sig"))
traceplot(output, pars=c("theta"))
traceplot(output, pars=c("ddb"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("beta_mu"))
stan_plot(output, pars=c("beta"))
stan_plot(output, pars=c("beta_sig"))
stan_plot(output, pars=c("theta"))
stan_plot(output, pars=c("ddb"))


############
print("comparing PA")
#pdf(paste("./plots/", model_name, "_comparison_PA.pdf", sep=""))
parameters = extract(output)
poseriorBeta1 = data.frame(parameters$beta[,,1])
poseriorBeta1 = reshape::melt(poseriorBeta1)
poseriorBeta2 = data.frame(parameters$beta[,,2])
poseriorBeta2 = reshape::melt(poseriorBeta2)
poseriorBeta3 = data.frame(parameters$beta[,,3])
poseriorBeta3 = reshape::melt(poseriorBeta3)
poseriorBeta4 = data.frame(parameters$beta[,,4])
poseriorBeta4 = reshape::melt(poseriorBeta4)



ggplot(poseriorBeta1, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("Beta 1") + theme(legend.position="none") + scale_x_continuous(limits=c(-15,30))
ggplot(poseriorBeta2, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("Beta 2") + theme(legend.position="none") + scale_x_continuous(limits=c(-15,30))
ggplot(poseriorBeta3, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("Beta 3") + theme(legend.position="none") + scale_x_continuous(limits=c(-15,30))
ggplot(poseriorBeta4, aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("Beta 4") + theme(legend.position="none") + scale_x_continuous(limits=c(-15,30))
ggplot(reshape::melt(data.frame(parameters$beta_mu)), aes(x=value, color=variable)) + geom_density() + geom_vline(aes(xintercept=mean(value))) + ggtitle("Beta Mu") + theme(legend.position="none") + scale_x_continuous(limits=c(-10,10))