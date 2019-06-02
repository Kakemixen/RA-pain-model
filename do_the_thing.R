# Programmed by Erling Ljunggren, cus that's important to know.

library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)

# set seed
set.seed(21111996) # such eggs of easter ^^

# source HDIofMCMC.R to calculate HDI
source("./HDIofMCMC.R")


### parameters for running
available_names = "
const_i
const_h
const_h_noT

lin_i
lin_h

log_i
log_h

exp_i

pow_i

"

model_name = "const_h_noT"

is_h = TRUE
has_ret = FALSE #doesnt affect atm
has_tau = FALSE

iters = 2000
warmups = 1000




# read the data file
dat_1 = read.table("./data/AllSubjectsProcessed.tsv", header=T, sep="\t")

allSubjs = unique(dat_1$SubjID)  # all subject IDs
N = length(allSubjs)                   # number of subjects
# T = table(dat_1)[1]     # number of trials per subject (=108)
# T = nrow(dat_1)     # number of trials per subject (=108)
T = 108
numIter = 108           # max number of iterations to find global minimum values
numPars = 3             # number of parameters


block_comment = "
So, of the basic I think we need
Trial - incrementign int
RiskType - 1|low, 2|med, 3|high
# TimeType - 1|early, 2|medium, 3|late - used as estimate for pain endured
RewardType - 1|low, 2|med, 3|high
Reward - float | actual reward
Shock - 0|no, 1|yes
ResponseType - 0|safe, 1|risky, NaN|no record
"

# data frames for fill-in
# -1 will be values for not used fields |  NaN
Tsubj = array(0, c(N))
RewardType = array(0, c(N,T))
RiskType = array(0, c(N,T))
ResponseType = array(0, c(N,T))
Shock = array(0, c(N,T))

#fill in with data
for (n in 1:N){
    subjdat = subset(dat_1, SubjID == allSubjs[n])
    trials = nrow(subjdat)
    Tsubj[n] = trials
    RiskType[n,1:trials] = subjdat$RiskType
    RewardType[n,1:trials] = subjdat$RewardType
    ResponseType[n,1:trials] = subjdat$ResponseType
    Shock[n,1:trials] = subjdat$Shock
}

dataList <- list(
    N       = N,
    T       = T, #108
    Tsubj   = Tsubj, # <= 108
    RiskType = RiskType,
    RewardType = RewardType,
    ResponseType =  ResponseType,
    Shock = Shock
)


# get pars vector
parameters = c("RiskAversion","PainAvoidance","log_lik","PredictedResponse")
if(has_tau){
      parameters = append(parameters, c("tau"))
}
if(is_h){
    parameters = append(parameters, c("mu_p", "sigma_p", "mu_RiskAversion","mu_PainAvoidance"))
    if(has_tau){
        parameters =append(parameters, c("mu_tau"))
    }
}
print("Parameters to fit")
print(parameters)

# run!
if (!file.exists(paste(model_name, "stanfit.rds",sep="_"))){
    print("fitting stan model")
    output = stan(paste("./models/", model_name, ".stan", sep=""),
          data = dataList,
          pars = parameters,
          iter = iters, warmup=warmups, chains=1, cores=1)
    saveRDS(output, paste(model_name, "stanfit.rds",sep="_"))
} else {
    print("recovering stan model")
    output = readRDS(paste(model_name, "stanfit.rds",sep="_"))
}


## traceplot
pdf(paste("./plots/", model_name, "_traceplot.pdf", sep=""))
traceplot(output, pars=c("RiskAversion"))
traceplot(output, pars=c("PainAvoidance"))
# traceplot(output, pars=c("tau"))
# traceplot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))


# posterior plots
pdf(paste("./plots/", model_name, "_posteriors.pdf", sep=""))
stan_plot(output, pars=c("RiskAversion"))
stan_plot(output, pars=c("PainAvoidance"))
# plot(density(extracted$PainAvoidance[,,3]), col="red")
# stan_dens(output, pars=c("tau"))
# stan_plot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))
# stan_dens(output, pars=c("mu_RiskAversion", "mu_PainAvoidance", "mu_tau"))
stan_plot(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))
stan_dens(output, pars=c("mu_RiskAversion", "mu_PainAvoidance"))


params = extract(output)
pdf(paste("./plots/", model_name, "_prediction.pdf", sep=""))

# Set up the vectors
true_gambling <- c("tSafe", "tGamble")
prediction <- c("pSafe","pGamble")

total_values <- c(0,0,0,0)

for(n in 1:N){
             #true: S, G  pred:
    pred_values <- c(0,0, # S
                     0,0) # G
    for(i in 1:Tsubj[n]){
        yPred <- params$PredictedResponse[iters-warmups,n,i]
        yTrue <- ResponseType[n,i]
        index = yPred*2 + yTrue + 1
        pred_values[index] = pred_values[index] + 1
    }
    total_values = total_values + pred_values

    # Create the data frame
    df <- expand.grid(true_gambling, prediction)
    df$value <- pred_values

    #Plot the Data
    g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
        theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste("Subject:",n)) +
        scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
    print(g)
}
#Plot the Data for all subjects together
df <- expand.grid(true_gambling, prediction)
df$value <- total_values
g <- ggplot(df, aes(Var1, Var2)) + geom_point(aes(size = value), colour = c("green", "red", "red", "green")) +
    theme(legend.position="none") + xlab("") + ylab("") + ggtitle(paste("All Subjects")) +
    scale_size_continuous(range=c(10,30)) + geom_text(aes(label = value))
print(g)

# pdf(paste("./plots/", model_name, "_likelihood.pdf", sep=""))
