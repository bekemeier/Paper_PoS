%% PostHoc Curves
% The script visualizes the time courses of post-hoc permutation-derived 
% t-statistics for Part-of-Speech (PoS) contrasts across cortical parcels. 
% It loads parcel-wise post-hoc results (verb vs. noun, adjective vs. noun, 
% verb vs. adjective) and reconstructs full brain × time matrices, 
% ensuring consistent parcel indexing across hemispheres and missing-data handling.

%For each parcel, the script extracts t-statistics and corresponding 
% p-values across time from the post-hoc permutation analyses and 
% reorganizes them into a standardized cortical format 
% (386 parcels × timepoints). These are then used to visualize temporal 
% dynamics of PoS-related effects at the parcel level.

%Finally, the script produces time-course plots for selected parcels, 
% showing the evolution of t-statistics over time for each PoS contrast. 

clear
thisDir = mfilename("fullpath")
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1)
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!

ft_defaults;
ft_hastoolbox('cellfunction',1)

dat_dir=fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir = fullfile(dat_dir,'t_stat_Ind_Pred'); 
savedir = fullfile(dat_dir,'Images'); 

filename = fullfile(basedirectory,'MOUSDATA','Data','Preproc','time.mat');
load(filename, 'time');
model_name = 'POS';

% load intact sentences
RA = zeros(numel(time),382);
RV = zeros(numel(time),382);
RVA = zeros(numel(time),382);
pval_A = zeros(numel(time),382);
pval_V = zeros(numel(time),382);
pval_VA = zeros(numel(time),382);
d = dir(fullfile(datadir,sprintf('PostHoc_Ind_Pred*%s.mat',model_name)));

numparc=length(d);
j=1;
for k = 1:382
    parcel_indx(k,1) = k;
    if  j>numparc
        RA(:,k) = zeros(numel(time),1);
        RV(:,k) = zeros(numel(time),1);
        RVA(:,k) = zeros(numel(time),1);
        pval_A(:,k) = 1;
        pval_V(:,k) = 1;
        pval_VA(:,k) = 1;
        
    elseif    parcel_indx(k,1) ~= str2num(d(j).name(26:28))
        RA(:,k) = zeros(numel(time),1);
        RV(:,k) = zeros(numel(time),1);
        RVA(:,k) = zeros(numel(time),1);
        pval_A(:,k) = 1;
        pval_V(:,k) = 1;
        pval_VA(:,k) = 1;
       
    else
       load(fullfile(datadir,d(j).name), 't_statA','t_statV','t_statVA','p_valsA','p_valsV','p_valsVA');
       RA(:,k) = t_statA(:);
       RV(:, k) = t_statV(:);
       RVA(:, k) = t_statVA(:);
       pval_A(:,k) = p_valsA(:);
       pval_V(:,k) = p_valsV(:);
       pval_VA(:,k) = p_valsVA(:);
       j=j+1;
      
    end
end
 
indx = 1:386;
indx([1 2 194 195]) = [];
% Create a field in the structure with the model name
adj = zeros(386, numel(time));
adj(indx, :) = RA';
verb = zeros(386, numel(time));
verb(indx, :) = RV';
verbadj = zeros(386, numel(time));
verbadj(indx, :) = RVA';
% do the same for the p-values
P_adj = zeros(386, numel(time));
P_adj(indx, :) = pval_A';
P_verb = zeros(386, numel(time));
P_verb(indx, :) = pval_V';
P_verbadj = zeros(386, numel(time));
P_verbadj(indx, :) = pval_VA';
%prep scatter plot for the sig effects
threshold=0.0167;



%%
parc_num=95; % for the LH (up to 191) add 2 to the original parcel number, from 192 always add 4 because of theis command: indx([1 2 194 195]) = [];
figure; 
plot(verb(parc_num,:),  'k','linewidth', 3);hold on
plot(adj(parc_num,:), '--k', 'linewidth', 2);hold on
plot(verbadj(parc_num,:), ':k', 'linewidth', 1.5);hold on
ylabel('t statistic')
xlim([0 109]);
set(gca, 'XTick', 1:12:length(time), 'XTickLabel', time(1:12:end))
xlabel('time (s)');
hold off
title(sprintf('Time course of t statistics at parcel %03d',parc_num-2))

% Define light brown as an RGB triplet
lightBrown = [0.71 0.40 0.11];  % adjust as needed

% Add dashed horizontal lines at y = 3.6 and y = -3.6 to mark the threshold
yline(3.6, '--', 'Color', lightBrown, 'LineWidth', 1);
yline(-3.6, '--', 'Color', lightBrown, 'LineWidth', 1);

legend({'Verb vs. Noun', 'Adjective vs. Noun', 'Verb vs. Adjective'},'Location','southoutside') 
ylim([-6 6]);
set(findobj(gcf,'type','axes'),'FontName','Arial','FontSize',12,'FontWeight','Bold', 'LineWidth', 2);

f=gcf;
f.WindowState='normal';
%%
plotname = fullfile(savedir, sprintf('%s',model_name), sprintf('PosHoc_Ind_Pred_Timecourse_%s_%03d_2',model_name,parc_num-2));
exportgraphics(f,sprintf('%s.pdf',plotname),'Resolution',300)
