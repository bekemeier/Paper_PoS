% This code writes the coefficient of determination (COD or R^2) for the
% original data and for the shuffled data from all scenarios into the same
% dataframe. The dimension of the dataframe: time X subject

% The script first loads the COD of the original data and averaged COD
% across 50 shuffled versions (dim: time x subject) from the model
% comparison scheme. For each parcel, the script concatenates the CODs from
% each scenario across the second dimension. The output of the script are
% two variables: CODall(109 x 100) and meanshufCODall (109 x 100).
% They are stored in the same file to be used later for the statistics.
clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);

%addpath  to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!
savedirectory = fullfile(basedirectory, 'MOUSDATA','Data');
ft_defaults;
ft_hastoolbox('cellfunction',1)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Beginning analysis script %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for p = 1:382
    parcel_indx = p;
    CODall = [];
    shufCODall = [];
    
    for scenario = 1:6
        dat_dir= fullfile(savedirectory,'Derivatives','Words_V_Sent');
        datadir = fullfile(dat_dir, sprintf('scenario%d', scenario), 'Stats');
        savedir = fullfile(dat_dir,'stats_allsce','Interactions');

        % Load COD and shuffled COD
        filename = fullfile(datadir, sprintf('mscca_sce%d_parcel%03d_Int_POS_Cont_FULLRED', scenario, parcel_indx));
        load(filename, 'COD','meanShufCOD');

         % Concatenate along subject dimension (2)
        CODall = [CODall, COD];                  % COD: [time × nSub]
        shufCODall = [shufCODall, meanShufCOD]; % shufCOD: [time × nSub]

        clear COD meanShufCOD;
    end

    % Save combined results
    filename = fullfile(savedir, sprintf('mscca_sceALL_parcel%03d_Int_POS_Cont_FULLRED', parcel_indx));
    save(filename, 'CODall', 'shufCODall', '-v7.3'); % v7.3 for large arrays
    clear CODall shufCODall
end