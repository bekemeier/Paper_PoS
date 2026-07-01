%% Ridge Regression
% The script loops through each scenario and parcel and fits subject-specific
% time-resolved encoding models using ridge regression. For each parcel, data
% from sentence and word-list conditions are combined and a full design matrix
% is constructed including ten predictors: (i) word length, (ii) lexical frequency,
% (iii) log probability, (iv) entropy, (v) ordinal position, (vi) context (−1/1),
% (vii) verb (0/1), (viii) adjective (0/1), and the interaction terms (ix) context×verb
% and (x) context×adjective.

% For each subject and timepoint, ridge regression is trained using five-fold
% cross-validation. The regularisation parameter (λ) is selected from a predefined
% grid by minimizing the mean squared error (MSE) across folds. The selected λ is
% then used to compute the final model error and coefficients.

% This procedure is repeated independently for each timepoint and subject,
% yielding time-resolved estimates of model performance and regression weights.

% In addition, a permutation-based null model is constructed by randomly
% permuting the rows of the design matrix 50 times within each subject and
% timepoint. For each permutation, ridge regression is refitted using the same
% cross-validation and λ-selection procedure. The resulting errors and coefficients
% are stored for each shuffle.

% All outputs are saved per parcel, including the original model estimates
% (err_full, err_red, coef_full, coef_red) and the corresponding permutation
% distributions (shuf_full, shuf_red, shufcoef_full, shufcoef_red).
clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);
datadir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent');
basedir =fullfile(basedirectory,'MOUSDATA','Data','Preproc');
%addpath  to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!
 
ft_defaults;
ft_hastoolbox('cellfunction',1)

%% ---- START PARALLEL POOL (ONCE) ----
maxNumCompThreads(1);
p = gcp('nocreate');
if isempty(p)
    parpool('local',4);   % fewer workers = more stability
end

%% ==================== SETTINGS ====================
nTime = 109;
nShuf = 50;
nCoef_full = 10;  % all+2 interactions
nCoef_red  = 8; % all individual predictors

lambdaList = [0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5];
%% ==================== SCENARIOS ====================
for scenario = 1:6

    load(fullfile(basedir,'subj_sce_info.mat'))
    sel  = strncmp(sce, num2str(scenario), 1);
    subj = subj(sel);
    subj = subj(strncmp(subj,'V',1));
    nSubj = numel(subj);
    clear sce

    fprintf('\nScenario %d — Huizeling interaction model\n', scenario)

    %% ==================== PARCEL LOOP (SERIAL) ====================
    for parcel_indx = 1:382

        fprintf('Scenario %d | Parcel %03d\n', scenario, parcel_indx)

        %% -------- Load & combine data --------
        fn = fullfile(datadir, sprintf('scenario%d',scenario), ...
            sprintf('mscca_sce%d_parcel%03d_tlck',scenario,parcel_indx));
        S = load(fn,'tlck');
        tlckSent = S.tlck;
        tlckSent.trialinfo.context = ones(height(tlckSent.trialinfo),1);

        fn = fullfile(datadir, sprintf('scenario%d',scenario), ...
            sprintf('mscca_sce%d_parcel%03d_tlckWL',scenario,parcel_indx));
        S = load(fn,'tlck');
        tlckWL = S.tlck;
        tlckWL.trialinfo.context = zeros(height(tlckWL.trialinfo),1);

        tlck = tlckSent;
        tlck.trial = cat(1, tlckSent.trial, tlckWL.trial);
        tlck.trialinfo = [tlckSent.trialinfo; tlckWL.trialinfo];

        clear tlckSent tlckWL S


     %% -------- Predictors (once per parcel) --------
        x1 = strlength(tlck.trialinfo.word);
        x2=9+log10(tlck.trialinfo.lexfreq);
        x2(~isfinite(x2)) = 0;
        x2 = (x2 - mean(x2)) ./ std(x2);

        x3 = -tlck.trialinfo.logprob ./ log(10);  % converts ln(p) to -log10(p)
        x3(~isfinite(x3)) = 0;                    % handle -Inf or NaN
        x3 = (x3 - mean(x3)) ./ std(x3);          % z-score


        x4=tlck.trialinfo.entropy;
        x4(~isfinite(x4)) = 0;
        x4 = (x4 - mean(x4)) ./ std(x4);

        x5=tlck.trialinfo.ordinal;
        x5 = (x5 - mean(x5)) ./ std(x5);
       % Context: effect coding (-1 = scrambled, +1 = intact)
        x6 = tlck.trialinfo.context;   % originally 0/1
        x6(x6 == 0) = -1;
        
        % POS dummy coding (reference = noun)
        x7 = tlck.trialinfo.Verb;      % 0/1
        x8 = tlck.trialinfo.Adj;       % 0/1
        
        % Interactions
        x9 = x6 .* x7;   % context × verb
        x10 = x6 .* x8;   % context × adjective
    
        X_full = [x1 x2 x3 x4 x5 x6 x7 x8 x9 x10];
        X_red  = [x1 x2 x3 x4 x5 x6 x7 x8];
        clear x1 x2 x3 x4 x5 x6 x7 x8 x9 x10

        %% -------- Allocate parcel results --------
        err_full  = zeros(nTime,nSubj);
        err_red   = zeros(nTime,nSubj);
        shuf_full = zeros(nTime,nShuf,nSubj);
        shuf_red  = zeros(nTime,nShuf,nSubj);

        coef_full = zeros(nTime,nCoef_full,nSubj);
        coef_red  = zeros(nTime,nCoef_red,nSubj);
        shufcoef_full = zeros(nTime,nCoef_full,nShuf,nSubj);
        shufcoef_red  = zeros(nTime,nCoef_red,nShuf,nSubj);

   %% ==================== SUBJECT PARFOR ====================
        parfor s = 1:nSubj
            rng(1000 + s, 'twister');
            

            err_full_s  = zeros(nTime,1);
            err_red_s   = zeros(nTime,1);
            shuf_full_s = zeros(nTime,nShuf);
            shuf_red_s  = zeros(nTime,nShuf);

            coef_full_s = zeros(nTime,nCoef_full);
            coef_red_s  = zeros(nTime,nCoef_red);
            shufcoef_full_s = zeros(nTime,nCoef_full,nShuf);
            shufcoef_red_s  = zeros(nTime,nCoef_red,nShuf);

            for tp = 1:nTime
                y = tlck.trial(:,s+3,tp);
                if any(isnan(y))
                    y = fillmissing(y,'constant',mean(y,'omitnan'));
                end

                n = numel(y);
                c = cvpartition(n,'KFold',5);

                [err_full_s(tp), coef_full_s(tp,:)] = ...
                    local_ridge_cv_coef(X_full,y,c,lambdaList);

                [err_red_s(tp), coef_red_s(tp,:)] = ...
                    local_ridge_cv_coef(X_red,y,c,lambdaList);

                for a = 1:nShuf
                    perm = randperm(n);
                    [shuf_full_s(tp,a), shufcoef_full_s(tp,:,a)] = ...
                        local_ridge_cv_coef(X_full(perm,:),y,c,lambdaList);
                    [shuf_red_s(tp,a), shufcoef_red_s(tp,:,a)] = ...
                        local_ridge_cv_coef(X_red(perm,:),y,c,lambdaList);
                end
            end

            err_full(:,s) = err_full_s;
            err_red(:,s)  = err_red_s;
            shuf_full(:,:,s) = shuf_full_s;
            shuf_red(:,:,s)  = shuf_red_s;

            coef_full(:,:,s) = coef_full_s;
            coef_red(:,:,s)  = coef_red_s;
            shufcoef_full(:,:,:,s) = shufcoef_full_s;
            shufcoef_red(:,:,:,s)  = shufcoef_red_s;
        end

        %% -------- Save parcel immediately --------
        outdir = fullfile(datadir, sprintf('scenario%d',scenario),'Models');
        if ~exist(outdir,'dir'); mkdir(outdir); end

        fnSave = fullfile(outdir, ...
            sprintf('mscca_sce%d_parcel%03d_Interaction_POSCont_FULLRED.mat', ...
            scenario, parcel_indx));

        save(fnSave, ...
            'err_full','err_red','shuf_full','shuf_red', ...
            'coef_full','coef_red','shufcoef_full','shufcoef_red', ...
            '-v7.3');

        clear tlck X_full X_red err_* coef_* shuf*
        if mod(parcel_indx,10)==0
            drawnow;
        end
    end
end

function [err, coefs] = local_ridge_cv_coef(X,y,c,lambdaList)
% Performs ridge regression with CV, returns MSE and average coefficients
nLambda = numel(lambdaList);
nFold = c.NumTestSets;
nCoef = size(X,2);
mse = zeros(nLambda,nFold);
coef_fold = zeros(nLambda,nCoef,nFold);

for i = 1:nLambda
    for f = 1:nFold
        tr = training(c,f);
        te = test(c,f);

        b = ridge(y(tr),X(tr,:),lambdaList(i),0);
        yhat = b(1) + X(te,:)*b(2:end);

        mse(i,f) = mean((yhat - y(te)).^2);
        coef_fold(i,:,f) = b(2:end);
    end
end

[~,idx] = min(mean(mse,2));
err = mean(mse(idx,:));
coefs = mean(coef_fold(idx,:,:),3);
end