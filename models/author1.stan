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

 

    // Group-level parameters.
    vector[n_pred]     beta_mu;   

}
model {

    // Model-generated data
    real    theta[n_obs];       // likelihood of take
    real    ddb[n_obs];         // distance from decision boundary

    // Hyperpriors (group-level priors)
    beta_mu ~ student_t(5, 0, 2); 
    
    // Likelihood
    for (n in 1:n_obs){
    
        // Generate theta/DDB for trial.
        theta[n] = inv_logit( dot_product( X[n], beta_mu ) );
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
    real log_lik[n_subj];
    
    for (i in 1:n_subj) {
      log_lik[i] = 0;
    }
    
    for (n in 1:n_obs){
        theta[n] = inv_logit( dot_product( X[n], beta_mu ) );
        ddb[n] = 0.25 - (theta[n]-0.5)^2;
        PointPosteriors[n] = exp( bernoulli_log( N[n], theta[n] ) );
        PredictedResponse[n] = bernoulli_rng(theta[n]);
        log_lik[ix[n]] += bernoulli_lpmf(N[n] | theta[n]);
    }

