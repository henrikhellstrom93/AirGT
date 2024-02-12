function f = qfunc_modified(gamma, Prob_np, n, P_max, sigma_z, sigma_h)
    % Returns the equation for q using P(n_p=i)
    % Input:
    % gamma = decision threshold (variable)
    % Prob_np = (n+1)-length vector that contains the 
    % probabilities P(n_p=i). Note that Prob_np is zero-indexed.
    % n = number of items
    % P_max = power constraint of user
    % sigma_z = additive noise ~ CN(0, 2*sigma_z^2)
    % sigma_h = channel ~ CN(0, 2*sigma_h^2)
    % Output:
    % f = equation of gamma that calculates q
    
    l = 1:n;
    % Factor in exponent of Eq. 22 in ICC paper
    gamma_coeff = -0.5 ./ (l * P_max * sigma_h^2 + sigma_z^2);
    % Factor in sum of Eq. 22 in ICC paper
    gamma_term = 1 - exp(gamma * gamma_coeff);
    % Eq. 22 in iCC paper
    f = max( exp(-gamma/(2*sigma_z^2)), gamma_term * Prob_np(2:n+1)/sum(Prob_np(2:n+1)) );
end

