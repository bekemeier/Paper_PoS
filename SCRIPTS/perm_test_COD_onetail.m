function results = perm_test_COD_onetail(COD, meanshufCOD, n_perm)

if nargin < 3
    n_perm = 5000;
end

[n_tp, n_sub] = size(COD);

%% -------------------------
%  Paired differences
% -------------------------
diffCOD = COD - meanshufCOD;   % (tp × sub)

%% -------------------------
%  Observed t-statistics
% -------------------------
Tobs = zeros(n_tp,1);

for tp = 1:n_tp
    x = diffCOD(tp,:);
    [~,~,~,stats] = ttest(x, 0, "Tail", "right");
    Tobs(tp) = stats.tstat;
end

%% -------------------------
%  Permutation null (sign-flipping)
% -------------------------
max_null = zeros(n_perm,1);

parfor k = 1:n_perm
    % Random ±1 sign for each subject
    signs = (rand(1,n_sub) > 0.5)*2 - 1;

    % Apply sign flip
    permData = diffCOD .* signs;

    % t-statistics across timepoints
    perm_t = zeros(n_tp,1);
    for tp = 1:n_tp
        x = permData(tp,:);
        if std(x) == 0
            perm_t(tp) = 0;
        else
            [~,~,~,stats] = ttest(x, 0, 'Tail','right');
            perm_t(tp) = stats.tstat;
        end
    end

    % max over time
    max_null(k) = max(perm_t);
end
%%  Confidence interval of max-null (95%)
% -------------------------
ci_null_max = prctile(max_null, 95);
%% -------------------------
%  Corrected p-values
% -------------------------
Pvals = zeros(n_tp,1);

for tp = 1:n_tp
    Pvals(tp) = mean(max_null >= Tobs(tp));
end

%% -------------------------
%  Output
% -------------------------
results.Tobs = Tobs;
results.Pvals = Pvals;
results.max_null = max_null;
results.ci_null_max = ci_null_max;
end
