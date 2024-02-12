function [pdf_np] = getPDFnp(n, d)
    % Computes the exact PDF P(n_p=i). See Eq. 21 in ICC paper.
    % Input:
    % n = number of devices
    % d = dimension of histogram
    % Output:
    % pdf_np = (n+1)-length vector that contains the probabilities P(n_p=i)

    i_list = 0:n;
    k_list = 0:d;
    pdf_np = zeros(length(i_list), 1);
    
    for i_idx = 1:length(i_list)
        i = i_list(i_idx)
    
        prob_participate = pdf('Binomial', i, n, k_list/d);
        prob_test_vector = pdf('Binomial', k_list, d, 0.5/n);
        prob_np_i = prob_participate*(prob_test_vector');
    
        pdf_np(i_idx) = prob_np_i;
    end
end

