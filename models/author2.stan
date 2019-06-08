data {
    // Metadata
    int n_obs;          // observations
    int n_pred;         // predictors
    int n_subj;         // subjects

    
    // Observations
    int N[n_obs];                 // choices
    int ix[n_obs];      // subject ix
    vector[n_pred] X[n_obs];      // predictors (risk, reward)
}
parameters {

    // Subject-level parameters.
    row_vector[n_pred] beta[n_subj];

    
    // Group-level parameters.
    real               beta_mu[n_pred];   
    real<lower=0>      beta_sig[n_pred];

}
model {

    // Model-generated data
    real    theta[n_obs];       // likelihood of take
    real    ddb[n_obs];         // distance from decision boundary



    beta_mu ~ student_t(5, 0, 2); 
    beta_sig ~ gamma(2,0.5);    // Equivalent to gamma(shape=2,scale=2) 
    
    
    for (i in 1:n_subj){
        for (j in 1:n_pred){
            beta[i,j] ~ student_t(5, beta_mu[j], beta_sig[j]); 
        }
    }
    
    // Likelihood
    for (n in 1:n_obs){
    
        // Generate theta/DDB for trial.
        theta[n] = inv_logit( dot_product( X[n], beta[ix[n]] ) );
        ddb[n] = 0.25 - (theta[n]-0.5)^2;
        
        // Model choice data.
        N[n] ~ bernoulli( theta[n] );
        
        
    }   
}
generated quantities {
    
    vector[n_obs] theta; 
    vector[n_obs] ddb; 
    vector[n_obs] PointPosteriors; // vector for computing log pointwise predictive density
    vector[n_obs] PredictedResponse;
    
    for (n in 1:n_obs){
        theta[n] = inv_logit( dot_product( X[n], beta[ix[n]] ) );
        ddb[n] = 0.25 - (theta[n]-0.5)^2;
        PointPosteriors[n] = exp( bernoulli_log( N[n], theta[n] ) );
        PredictedResponse[n] = bernoulli_rng(theta[n]);
    }
}