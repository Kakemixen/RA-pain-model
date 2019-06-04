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
    // hyper raw
    vector[6] mu_p;
    vector<lower=0>[6] sigma_p;

    //individ raw
    vector[N] RiskAversion_pr;
    matrix[N, 3] PainAvoidance_pr;
    vector[N] PainRetention_pr;
    vector[N] tau_pr;
}
transformed parameters {
    vector<lower=0, upper=5>[N] RiskAversion;
    matrix<lower=0, upper=5>[N, 3] PainAvoidance;
    vector<lower=0, upper=5>[N] PainRetention;
    vector<lower=0, upper=20>[N] tau;

    for(i in 1:N){
        RiskAversion[i] = Phi_approx( mu_p[1] + sigma_p[1] * RiskAversion_pr[i] ) * 5;
        //Pain avoidance is matrix
        for (j in 1:3){
            PainAvoidance[i,j] = Phi_approx( mu_p[j+1] + sigma_p[j+1] * PainAvoidance_pr[i,j] ) * 5;
        }
        PainRetention[i] = Phi_approx( mu_p[5] + sigma_p[5] * PainRetention_pr[i] ) * 5;
        tau[i] = Phi_approx( mu_p[6] + sigma_p[6] * tau_pr[i] ) * 20;
    }
}
model {
    // priors
    // hyper
    mu_p ~ normal(0,1);
    sigma_p ~ normal(0,1);
    // Individual
    RiskAversion_pr    ~ normal(0, 1);
    PainRetention_pr   ~ normal(0, 1);
    to_vector(PainAvoidance_pr) ~ normal(0, 1);
    tau_pr    ~ normal(0, 1);

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

            evGamble = pow(RewardType[i,t]*0.33, RiskAversion[i]) - PainAvoidance[i,RiskType[i,t]] * exp(n_shocks * PainRetention[i]);

            pGamble  = inv_logit(tau[i] * (evGamble - evSafe));

            ResponseType[i,t] ~ bernoulli(pGamble);

            // update shocks, RL?
            if(ResponseType[i,t] == 1)
                n_shocks += Shock[i,t];
        }
    }
}

generated quantities {
    real mu_RiskAversion;
    vector[3] mu_PainAvoidance;
    real mu_PainRetention;
    real mu_tau;

    // For posterior predictive check
    real PredictedResponse[N, T];
    real log_lik[N];

    mu_RiskAversion = Phi_approx(mu_p[1]) * 5;
    mu_PainAvoidance = Phi_approx(mu_p[2:4]) * 5;
    mu_PainRetention = Phi_approx(mu_p[5]) * 5;
    mu_tau = Phi_approx(mu_p[6]) * 20;



    // Set all posterior predictions to 0 (avoids NULL values)
    for (i in 1:N) {
        for (t in 1:T) {
            PredictedResponse[i, t] = -1;
        }
    }

    { // local section, this saves time and space
        for (i in 1:N) {
            int n_shocks = 0; //remember to reset for each sequence of trial / subject
            log_lik[i] = 0;
            for (t in 1:Tsubj[i]) {
                real evSafe;
                real evGamble;
                real pGamble;

                evSafe     = pow(0.01, RiskAversion[i]);
                evGamble   = pow(RewardType[i,t]*0.33, RiskAversion[i]) - PainAvoidance[i,RiskType[i,t]] * exp(n_shocks * PainRetention[i]);
                pGamble    = inv_logit(tau[i] * (evGamble - evSafe));
                log_lik[i] += bernoulli_lpmf(ResponseType[i, t] | pGamble);

                // generate posterior prediction for current trial
                PredictedResponse[i, t] = bernoulli_rng(pGamble);
                // update shocks, RL?
                if(ResponseType[i,t] == 1)
                    n_shocks += Shock[i,t];
            }
        }
    }
}
