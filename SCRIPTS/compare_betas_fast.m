function [t_stat, p_vals, sig_times] = compare_betas_fast(betas1, betas2, n_perm, alpha, time)

% betas1, betas2: [timepoints × subjects]
% n_perm: number of permutations
% alpha: significance threshold
% time (optional): physical time vector

if nargin < 5
    time = 1:size(betas1,1);
end

[nTimepoints, nSubjects] = size(betas1);

%% ----------------------------------------
% Observed paired t-statistics (vectorized)
% ----------------------------------------
diffDATA = betas1 - betas2;              % [tp × subj]
mu  = mean(diffDATA, 2);                 % [tp × 1]
sd  = std(diffDATA, 0, 2);               % [tp × 1]
t_stat = (mu ./ (sd ./ sqrt(nSubjects)))';  % row vector

%% ----------------------------------------
% Permutation null (sign flipping)
% ----------------------------------------
max_t_null = zeros(n_perm,1);

parfor perm = 1:n_perm
    flipSigns = (rand(1,nSubjects) > 0.5)*2 - 1;
    permData  = diffDATA .* flipSigns;   % implicit expansion

    mu_p = mean(permData, 2);
    sd_p = std(permData, 0, 2);
    t_p  = mu_p ./ (sd_p ./ sqrt(nSubjects));

    max_t_null(perm) = max(abs(t_p));
end

%% ----------------------------------------
% Corrected p-values
% ----------------------------------------
p_vals = mean(max_t_null >= abs(t_stat), 1);

%% ----------------------------------------
% Significant timepoints
% ----------------------------------------
sig_times = time(p_vals < alpha);

end