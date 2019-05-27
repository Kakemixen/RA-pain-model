/*
parameters:
pain avoidance (dependent on total previous pain)
subjective value
tau, same as ra | set as 1 maybe?
*/

data {
    int<lower=1> T;

    // Risk Aversion Part
    int<lower=1, upper=3> RiskType[T];
    int<lower=1, upper=3> RewardType[T];
    int<lower=0, upper=1> ResponseType[T];

    // Reinforcement Learning part
    // real<lower=0> Reward[T];
    int<lower=0, upper=1> Shock[T];
}
transformed data {
}
parameters {
    real<lower=0, upper=2> RiskAversion;
    real<lower=0, upper=2> PainAvoidance;
    real<lower=0, upper=10> tau;
}
transformed parameters {
}
model {
    // counting shocks
    int n_shocks = 0; //remember to reset for each sequence of trial / subject
    //
    RiskAversion    ~ uniform(0, 2);
    PainAvoidance ~ uniform(0, 2);
    tau    ~ uniform(0, 10);


    for (t in 1:T) {

        real evSafe;
        real evGamble;
        real pGamble;

        evSafe   = pow(0.01, RiskAversion);

        if (RiskType[t] == 1)
            evGamble = 0.9 * (pow(RewardType[t], RiskAversion)) + 0.1 * (pow(RewardType[t], RiskAversion) - PainAvoidance * n_shocks);
        else if (RiskType[t] == 2)
            evGamble = 0.5 * (pow(RewardType[t], RiskAversion)) + 0.5 * (pow(RewardType[t], RiskAversion) - PainAvoidance * n_shocks);
        else
            evGamble = 0.1 * (pow(RewardType[t], RiskAversion)) + 0.9 * (pow(RewardType[t], RiskAversion) - PainAvoidance * n_shocks);

        pGamble  = inv_logit(tau * (evGamble - evSafe));

        ResponseType[t] ~ bernoulli(pGamble);

        // update shocks
        if(ResponseType[t] == 1)
            n_shocks += Shock[t];
    }
}
/*
generated quantities{
    // For posterior predictive check
    real y_pred[T];

    for (t in 1:T) {
        real evSafe;
        real evGamble;
        real pGamble;

        evSafe     = pow(cert[t], RiskAversion);
        evGamble   = 0.5 * (pow(gain[t], RiskAversion) + PainAvoidance * pow(loss[t], RiskAversion));
        pGamble    = inv_logit(tau * (evGamble - evSafe));

        // generate posterior prediction for current trial
        y_pred[t] = bernoulli_rng(pGamble);
    }
}
*/

