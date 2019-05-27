# Programmed by Erling Ljunggren, cus that's important to know.

library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)

# set seed
set.seed(21111996) # such eggs of easter ^^

# source HDIofMCMC.R to calculate HDI
source("./HDIofMCMC.R")

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
Tsubj = array(0, c(N))
RewardType = array(0, c(N,T))
RiskType = array(0, c(N,T))
ResponseType = array(0, c(N,T))
Shock = array(0, c(N,T))

print("subset")
print(length(subset(dat_1, SubjID == 1)[,1]))
print("riskty")
print(length(RiskType[1,]))
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

    #  # Risk Aversion Part
    # int<lower=1, upper=3> RiskType[T];
    # int<lower=1, upper=3> RewardType[T];
    # int<lower=0, upper=1> ResponseType[T];
    #
    # # Reinforcement Learning part
    # # real<lower=0> Reward[T];
    # int<lower=0, upper=1> Shock[T];
    #
dataList <- list(
    N       = N,
    T       = T, #108
    Tsubj   = Tsubj,
    RiskType = RiskType,
    RewardType = RewardType,
    ResponseType =  ResponseType,
    Shock = Shock
)



# output = stan("./model_1.stan",
#           data = dataList,
#           pars = c("RiskAversion", "PainAvoidance", "tau"),
#           iter = 4000, warmup=2000, chains=4, cores=4)
# # run!
# if (!file.exists("stanfit.rds")){
#     print("fitting stan model")
#     output = stan("./model_1.stan",
#               data = dataList,
#               pars = c("RiskAversion", "PainAcceptance", "tau"),
#               iter = 2000, warmup=1000, chains=4, cores=4)
#     saveRDS(output, "stanfit.rds")
# } else {
#     print("recovering stan model")
#     output = readRDS("stanfit.rds")
# }
#
# pdf("traceplot.pdf")
# traceplot(output)
# pdf("posteriors.pdf")
# stan_dens(output)
#
