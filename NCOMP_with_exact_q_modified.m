clear
clc

%% parameter setup
d = 1e4;                % dimension of each item
n = 10;                 % number of items
alpha = 0.5;            % parameter for p
p = alpha/n;            % group selection parameter
delta = 1/4;            % parameter for histogram estimation error probability
p_error = d^(-delta);   % histogram estimation error probability bound

P_max = 1;              % power constraint of user
sigma_h = 1;            % channel ~ CN(0, 2*sigma_h^2)
SNR_dB = 25;             % SNR = P_max*sigma_h^2 / sigma_z^2
sigma_z = sqrt(P_max * sigma_h^2 / 10^(SNR_dB/10));    % additive noise ~ CN(0, 2*sigma_z^2)
T = 500;     % # group tests

%% calculate q
tic

% Calculate pdf of n_p
Prob_np = getPDFnp(n, d);   % Prob_np(i) = P(n_p=i-1)

% set equation for q
gamma_upper_bound = ...
    2*(5*P_max * sigma_h^2 + sigma_z^2); %Just loose bound to init fminbnd
% Eq. 23 in ICC paper
gamma = fminbnd(@(x)qfunc_modified(x, ...
    Prob_np, n, P_max, sigma_z, sigma_h), 0, gamma_upper_bound);
% Select q based on best gamma (found through fminbnd)
q = qfunc_modified(gamma, Prob_np, n, P_max, sigma_z, sigma_h)

time = toc;
disp(['It took ', num2str(time, 2), ' seconds to calculate q.'])

%% setup NCOMP
cap_delta = sqrt(delta) * exp(-alpha) * (1-2*q) ...
    / (sqrt(delta) + sqrt(delta+1)) / q; % See Eq. (4) in ICC paper
beta = 2*exp(1)*(sqrt(delta)+sqrt(delta+1))^2*log(2)/(1-exp(-2))/(1-2*q)^2;
T_uppbound = beta*n*log(d)/log(2); % See Eq. (5) in ICC paper


num_mc = 100;

num_errors = 0;

%% run NCOMP

tic
for i = 1:num_mc
    i
    % Generate items
    x = zeros(d,1);
    for j = 1:n
        ind = randi(d);
        x(ind) = x(ind)+1;
    end
    
    % Generate group-testing matrix (Bernoulli method)
    M = false(T, d);
    for k = 1:T
        M(k, :) = logical(binornd(1, p, [1,d]));
    end
    
    % Get noiseless test results
    n_p = M*x;

    % Calculate y with random channel and noise
    y = zeros(T, 1);
    for j = 1:T
        h = sigma_h*(randn(1, n_p(j)) + 1j*randn(1, n_p(j)));
        z = sigma_z*(randn(1) + 1j*randn(1));
        y(j) = sum(P_max^0.5 * h) + z;
    end
    
    % Solve decision problem with threshold gamma
    y_dec = abs(y).^2>=gamma;

    % Find the vector lengths ||m_j||_1.
    % The variable "indicator" is a d-length vector containing the
    % norms of ||m_j||_1 for all j=1:d.
    indicator = sum(M, 1);
    
    % Find the size of the matching set.
    % The 'matching' variable is a d-length vector containing the
    % cardinalities of S_j for all j=1:d.
    M_and_y = bitand(M, y_dec);
    matching = sum(M_and_y, 1);

    % Estimate x (Eq. 3 in ICC paper)
    x_hat = ( matching >= indicator*(1-q*(1+cap_delta)) )';
    
    % Results
    idx_cmp = xor(x, x_hat);

    if sum(idx_cmp) ~= 0
        num_errors = num_errors + 1;
    end
end
time = toc;

frac_errors = num_errors/num_mc;

disp(['It took ', num2str(time, 2), ' seconds to complete the simulations.'])
disp([num2str(frac_errors*100, 3), '% of the monte-carlo trials had an error.'])