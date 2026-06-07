function aic  = aic_data(simulated_data, observed_data, k)

% Step 1: Fit a multivariate normal distribution to the simulated data
mean_estimated = mean(simulated_data);
cov_estimated = cov(simulated_data);

% Step 2: Calculate the log-likelihood for the observed data
log_likelihood = sum(log(mvnpdf(observed_data, mean_estimated, cov_estimated)));

aic = 2*k-2*log_likelihood;
% bic = k*log(length(observed_data))-2*log_likelihood;
