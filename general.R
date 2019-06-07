# Programmed by Erling Ljunggren, cus that's important to know.

library(rstan)
library(loo)
library(ggplot2)
rstan_options(auto_write = TRUE)

# set seed
set.seed(21111996) # such eggs of easter ^^

sample_model <- function(model_name, dataList, paramList, iterations, warmups, chains){
    # run!
    if (!file.exists(paste("./stanfits/", model_name, ".rds",sep=""))){
        print("fitting stan model")
        output = stan(paste("./models/", model_name, ".stan", sep=""),
              data = dataList,
              pars = paramList,
              iter = iterations, warmup=warmups, chains=chains, cores=chains)
        saveRDS(output, paste("./stanfits/", model_name, ".rds",sep=""))
    } else {
        print("recovering stan model")
        output = readRDS(paste("./stanfits/", model_name, ".rds",sep=""))
    }
    return(output)
}

get_dataList <- function(){
    # read the data file
    dat_1 = read.table("./data/AllSubjectsProcessed.tsv", header=T, sep="\t")

    allSubjs = unique(dat_1$SubjID)  # all subject IDs
    N = length(allSubjs)                   # number of subjects
    # T = table(dat_1)[1]     # number of trials per subject (=108)
    # T = nrow(dat_1)     # number of trials per subject (=108)
    T = 108
    numIter = 108           # max number of iterations to find global minimum values
    numPars = 3             # number of parameters


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
    return(dataList)
}

getmode <- function(v) {
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

PPC <- function(output, dataList){
    print("running PPC")
    pdf(paste("./plots/", model_name, "_prediction.pdf", sep=""))
    params = rstan::extract(output)

    # Set up the vectors
    true_gambling <- c("tSafe", "tGamble")
    prediction <- c("pSafe","pGamble")

    total_values <- c(0,0,0,0)

    for(n in 1:dataList$N){
                 #true: S, G  pred:
        pred_values <- c(0,0, # S
                         0,0) # G
        for(i in 1:dataList$Tsubj[n]){
            yPred <- getmode(params$PredictedResponse[,n,i])
            yTrue <- dataList$ResponseType[n,i]
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
}

getBIC <- function(ll, num_params, samples){
   (-2) * mean(ll) + log(samples) * num_params
}

BIC <- function(output, dataList, num_params, samples){
    formula = "ln(n)k -2ln(L)"
    print("running BIC")
    pdf(paste("./plots/", model_name, "_BIC.pdf", sep=""))
    params = rstan::extract(output)
    individ_BIC = rep(0, dataList$N)
    for(n in 1:dataList$N){
        individ_BIC[n] = getBIC(params$log_lik, num_params, samples)
    }
    df_individ_BIC = data.frame(BIC = individ_BIC, id = 1:dataList$N)
    g <- ggplot(data = df_individ_BIC, mapping = aes(x=id, y=BIC)) + geom_point() + ggtitle(paste("BIC per subject | total BIC:", sum(individ_BIC)))
    print(g)
}

LOOIC <- function(output){
    # Extract pointwise log-likelihood and compute LOO
    log_lik <- extract_log_lik(output, merge_chains = FALSE)

    # as of loo v2.0.0 we can optionally provide relative effective sample sizes
    # when calling loo, which allows for better estimates of the PSIS effective
    # sample sizes and Monte Carlo error
    r_eff <- relative_eff(exp(log_lik))

    looic <- loo(log_lik, r_eff = r_eff, cores = 2)
    print(looic)
}
