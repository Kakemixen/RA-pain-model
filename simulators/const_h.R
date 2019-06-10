library(boot)

# Simulation parameters
seed <- 211196
num_subjs  <- 36 # number of subjects | just setting a number larger than subjects in trials

# Set seed
set.seed(seed)   # always set a seed number for this homework!

# True parameters
simul_pars <- data.frame(
                         RiskAversion = abs(rnorm(num_subjs, 0.50, 0.8)),
                         PainAvoidance_low = abs(rnorm(num_subjs, 0.10, 0.20)),
                         PainAvoidance_med = abs(rnorm(num_subjs, 0.60, 0.30)),
                         PainAvoidance_high = abs(rnorm(num_subjs, 1.20, 0.40)),
                         tau = abs(rnorm(num_subjs, 10.00, 4.0)))


# read the data file | because we'll just use the presented conditions cus they are equally devided
dat_1 = read.table("../data/AllSubjectsProcessed.tsv", header=T, sep="\t")

# For storing simulated choice data for all subjects
all_data = data.frame(Trial = dat_1$Trial,
                     RiskType = dat_1$RiskType,
                     RewardType = dat_1$RewardType,
                     ResponseType = dat_1$ResponseType,
                     Reward = dat_1$Reward,
                     Shock = dat_1$Shock,
                     SubjID = dat_1$SubjID)

responses = rep(0, nrow(dat_1))
for (i in 1:nrow(dat_1)) {
    # Individual-level (i.e. per subject) parameter values
    RiskAversion <- simul_pars$RiskAversion[dat_1$SubjID[i]]
    PainAvoidance <- c(simul_pars$PainAvoidance_low[dat_1$SubjID[i]], simul_pars$PainAvoidance_med[dat_1$SubjID[i]], simul_pars$PainAvoidance_high[dat_1$SubjID[i]])
    tau <- simul_pars$tau[dat_1$SubjID[i]]

    evSafe  <- 0.01^RiskAversion

    evGamble <- (all_data$RewardType[i]*0.33)^RiskAversion - log(PainAvoidance[all_data$RiskType[i]] + 1)

    pGamble <- inv.logit(tau * (evGamble - evSafe))

    ResponseType <- rbinom(1,1,pGamble)

    # Append current subject with all subjects' data
    responses[i] = ResponseType
}

all_data$ResponseType = responses

# Write out data
write.table(all_data, file = "simul_data_const_h.txt", row.names = F, col.names = T, sep = "\t")
# write out parameters
write.table(simul_pars[-c(14),], file ="simul_param_const_h.txt", row.names=F, col.names = T, sep = "\t")
write.table(
            data.frame(
              RiskAversion = 0.5,
              PainAvoidance_low = 0.1,
              PainAvoidance_med = 0.5,
              PainAvoidance_high = 0.9,
              tau = 10
            ),
            file ="simul_param_hyp_const_h.txt", row.names=F, col.names = T, sep = "\t")
