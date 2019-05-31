/*
parameters:
pain avoidance (dependent on total previous pain)
subjective value
tau, same as ra | set as 1 maybe?
*/

data {
    int<lower=1> T;
    int<lower=1> N;
    int<lower=1> Tsubj[N];

    // Risk Aversion Part
    int<lower=0, upper=3> RiskType[N,T];
    int<lower=0, upper=3> RewardType[N,T];
    int<lower=0, upper=1> ResponseType[N,T];

    // Reinforcement Learning part
    // real<lower=0> Reward[T];
    int<lower=0, upper=1> Shock[N,T];
}
transformed data {
}
parameters {
    vector<lower=0, upper=3>[N] RiskAversion;
    vector<lower=0, upper=5>[N] PainAvoidance;
    vector<lower=0, upper=10>[N] tau;
}
transformed parameters {
}
model {
    // priors
    RiskAversion    ~ uniform(0, 3);
    PainAvoidance ~ uniform(0, 5);
    tau    ~ uniform(0, 10);

    // model calculation
    for (i in 1:N){
    // counting shocks
    int n_shocks = 0; //remember to reset for each sequence of trial / subject

        for (t in 1:Tsubj[i]) {

            // Risk Aversion
            real evSafe;
            real evGamble;
            real pGamble;

            evSafe   = pow(0.01, RiskAversion[i]);

            if (RiskType[i,t] == 1)
                evGamble = pow(RewardType[i,t]*0.33, RiskAversion[i]) - 0.1 * log( PainAvoidance[i] * n_shocks + 1);
            else if (RiskType[i,t] == 2)
                evGamble = pow(RewardType[i,t]*0.33, RiskAversion[i]) - 0.5 * log( PainAvoidance[i] * n_shocks + 1);
            else
                evGamble = pow(RewardType[i,t]*0.33, RiskAversion[i]) - 0.9 * log( PainAvoidance[i] * n_shocks + 1);

            pGamble  = inv_logit(tau[i] * (evGamble - evSafe));

            ResponseType[i,t] ~ bernoulli(pGamble);

            // update shocks, RL?
            if(ResponseType[i,t] == 1)
                n_shocks += Shock[i,t];
        }
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

