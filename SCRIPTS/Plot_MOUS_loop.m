%% Onset Mapping of Significant Interaction Effects
% This script computes and visualizes the temporal onset of significant
% interaction effects across cortical parcels. For each model, parcel-level
% t-statistics are loaded together with a list of FDR-corrected significant parcels.
% T-values are thresholded using the significance threshold
% (t_thresh) and restricted to parcels surviving FDR correction. The resulting
% time-resolved binary significance matrix is used to identify, for each parcel,
% the first timepoint at which a significant effect emerges after stimulus onset
% (t >= 0). Onset latency is defined as the earliest significant sample and is
% converted to milliseconds.
% Cortical surface geometry is loaded from the Conte69 atlas and inflated
% surface reconstruction. Left and right hemispheres are extracted, and
% parcel-level onset values are projected onto vertices according to the
% parcellation scheme.
% Separate surface meshes are constructed for each hemisphere and visualized
% using vertex-based color mapping. The resulting figures display spatial
% distributions of onset latencies for significant interaction effects across
% lateral and medial views of both hemispheres. A shared colorbar encodes
% onset time in milliseconds.

clear
thisDir = mfilename("fullpath");
basedirectory = thisDir(1:regexp(thisDir,'MOUSDATA')-1);
datadir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent');
basedir =fullfile(basedirectory,'MOUSDATA','Data','Preproc');
%addpath  to the fieldtrip folder
addpath /mnt/smbdir/fieldtrip-20230402; %check!!!
 
ft_defaults;
ft_hastoolbox('cellfunction',1)clear

data_dir = fullfile(basedirectory,'MOUSDATA','Data','Derivatives','Words_V_Sent','stats_allsce');
datadir = fullfile(data_dir,'t_stat_Interactions'); %interactions
savedir = fullfile(data_dir,'Images_Interactions'); % interactions
% % Load necessary data
filename = fullfile(basedir,'MOUSDATA','Data','Preproc','time.mat');
load(filename, 'time');
t0_idx = find(time >= 0);
nTime  = numel(time);
modelname = {'Cont'};%};'Entropy', 'Index','Lexfreq','Word'
%%
for j=1:numel(modelname)
     model_name = modelname{j};
disp(['Processing: ', model_name]);


% Initialize results matrix
R = zeros(nTime,382);


filename=fullfile(data_dir, sprintf('Interaction_FDR_%s_FMperm.mat',model_name)); % IndPred for individual predictors, Interaction
load(filename, 'all_t','significant_parcels','t_thresh')


% Set threshold for significant t-values
sig_threshold = t_thresh;
%% Load t-values into a matrix
R(t0_idx, :) = all_t;
% Apply significance threshold
R(R < sig_threshold) = 0;
R(:, ~ismember(1:size(R, 2), significant_parcels)) = 0;
%%
indx = 1:386;
indx([1 2 194 195]) = [];
dat_onset1 = zeros(386, numel(time));
 
dat_onset1(indx,:) = R';

%% ------------------------------------------------------------------------
%            Compute ONSET of significant effect per parcel
% -------------------------------------------------------------------------

onset_map = nan(386,1);   % onset in ms; NaN = no effect


for k = 1:386
   
    sig_idx = find(dat_onset1(k, t0_idx) ~= 0);

    if ~isempty(sig_idx)
        onset_idx = t0_idx(sig_idx(1));   % first significant sample ≥ 0
        onset_map(k) = time(onset_idx) * 1000;  % s → ms
    end
end

%% --- FIXED COLOR RANGE: 0–800 ms ---
onset_min = 0;
onset_max = 800;

valid = isfinite(onset_map);
onset_map(valid & onset_map < onset_min) = onset_min;
onset_map(valid & onset_map > onset_max) = onset_max;

%% Create a fixed colormap (256 steps)
cmap_onset = parula(256);  % early = dark blue, late = yellow/orange

% Normalize onset values into [1,256]
norm_onset = nan(size(onset_map));

norm_onset(valid) = round( ...
    1 + (onset_map(valid) - onset_min) ./ ...
    (onset_max - onset_min) * 255 );

norm_onset(~valid) = 1;     % dummy value, will never be used
norm_onset(norm_onset < 1)   = 1;
norm_onset(norm_onset > 256) = 256;
%% Load brain atlas & adjust positioning
load atlas_conte69_8196reg_LR_brodmann_subparc.mat
load cortex_inflated_8196reg.mat
atlas.pos = sourcemodel.pnt;
atlas.pos(4099:end,:) = atlas.pos(4099:end,:) * diag([-1 -1 1]);
m1 = mean(atlas.pos(1:4098,:)); atlas.pos(1:4098,:) = atlas.pos(1:4098,:) - m1;
m2 = mean(atlas.pos(4099:end,:)); atlas.pos(4099:end,:) = atlas.pos(4099:end,:) - m2;
atlas.pos(1:4098,2) = atlas.pos(1:4098,2) + 110;
atlas.pos(4099:end,2) = atlas.pos(4099:end,2) - 110;

left_hemi_indices = (atlas.brainstructure == 1);  % Assuming 1 indicates left hemisphere

% Filter vertices and faces for the left hemisphere
left_vertices = atlas.pos(left_hemi_indices, :);

% Find the faces that are entirely within the left hemisphere
left_faces_mask = all(left_hemi_indices(atlas.tri), 2);
left_faces = atlas.tri(left_faces_mask, :);

% Adjust face indices to correspond to the subset of vertices
[~, ~, new_face_indices] = unique(left_faces);
left_faces = reshape(new_face_indices, size(left_faces));
left_parcels = atlas.parcellation(1:4098);

left_parcellabels = atlas.parcellationlabel(1:193);

% Create a mesh structure for the left hemisphere
mesh_left_hemi.pos = left_vertices;
mesh_left_hemi.tri = left_faces; 
mesh_left_hemi.unit = atlas.unit;
mesh_left_hemi.parcellation = left_parcels;
mesh_left_hemi.parcellationlabel = left_parcellabels;
mesh_left_hemi.brainstructure =atlas.brainstructure(1:4098,:);
mesh_left_hemi.brainstructurelabel = atlas.brainstructurelabel{1,1};
mesh_left_hemi.inside = atlas.inside(1:4098,:);

% Determine the indices for the right hemisphere
right_hemi_indices = (atlas.brainstructure == 2);  % Assuming 1 indicates left hemisphere

% Filter vertices and faces for the left hemisphere
right_vertices = atlas.pos(right_hemi_indices, :);

% Find the faces that are entirely within the left hemisphere
right_faces_mask = all(right_hemi_indices(atlas.tri), 2);
right_faces = atlas.tri(right_faces_mask, :);

% Adjust face indices to correspond to the subset of vertices
[~, ~, right_face_indices] = unique(right_faces);
right_faces = reshape(right_face_indices, size(right_faces));
right_parcels = atlas.parcellation(4099:end);

right_parcellabels = atlas.parcellationlabel(194:end);

% Create a mesh structure for the left hemisphere
mesh_right_hemi.pos = right_vertices;
mesh_right_hemi.tri = right_faces; 
mesh_right_hemi.unit = atlas.unit;
mesh_right_hemi.parcellation = right_parcels;
mesh_right_hemi.parcellationlabel = right_parcellabels;
mesh_right_hemi.brainstructure =atlas.brainstructure(4099:end,:);
mesh_right_hemi.brainstructurelabel = atlas.brainstructurelabel{1,2};
mesh_right_hemi.inside = atlas.inside(4099:end,:);


% Get number of vertices in each hemisphere
num_vertices_LH = length(left_parcels);
num_vertices_RH = length(right_parcels);


% Initialize matrices for the brain vertices and transparency
vertex_colors_LH = zeros(num_vertices_LH, 3); % RGB colors for the left hemisphere
vertex_colors_RH = zeros(num_vertices_RH, 3); % RGB colors for the right hemisphere
vertex_alpha_LH = zeros(num_vertices_LH, 1);  % Transparency for the left hemisphere
vertex_alpha_RH = zeros(num_vertices_RH, 1);  % Transparency for the right hemisphere

%% --- MAP COLORS TO LEFT HEMISPHERE ---
for i = 1:num_vertices_LH
    parcel_idx = mesh_left_hemi.parcellation(i);

    if parcel_idx > 0 && parcel_idx <= 386
        idx_color = norm_onset(parcel_idx);
        if onset_map(parcel_idx) > 0
            vertex_colors_LH(i,:) = cmap_onset(idx_color, :);
        else
            vertex_colors_LH(i,:) = [1 1 1]; % white for no effect
        end
    else
        vertex_colors_LH(i,:) = [1 1 1];
    end
end

%% --- MAP COLORS TO RIGHT HEMISPHERE ---
for i = 1:num_vertices_RH
    parcel_idx = mesh_right_hemi.parcellation(i);

    if parcel_idx > 0 && parcel_idx <= 386
        idx_color = norm_onset(parcel_idx);
        if onset_map(parcel_idx) > 0
            vertex_colors_RH(i,:) = cmap_onset(idx_color, :);
        else
            vertex_colors_RH(i,:) = [1 1 1]; % white
        end
    else
        vertex_colors_RH(i,:) = [1 1 1];
    end
end

%% Plotting the results 
figure('Name', 'Significant Effects Left Lateral');
ft_plot_mesh(mesh_left_hemi, 'vertexcolor', vertex_colors_LH, 'edgecolor', 'none');
view([-90 0]);
lighting gouraud;
material dull;
h1 = light('position',[-100 0 50]); 
h2 = light('position',[-100 0 -50]);
set(gcf,'color','w');
f=gcf;
f.WindowState='maximized';
% Add the first colorbar for "Onset of the effect"
colormap(cmap_onset);
c = colorbar;
c.Label.String = 'Onset (ms)';
c.Ticks = linspace(0,1,9);
c.TickLabels = 0:100:800;
c.FontSize = 12;
plotname = fullfile(savedir, sprintf('%s',model_name), sprintf('onset%s_brain_LateralLeft', model_name));
exportgraphics(f,sprintf('%s.png',plotname),'Resolution',300)

figure('Name', 'Significant Effects Right Lateral')
ft_plot_mesh(mesh_right_hemi, 'vertexcolor', vertex_colors_RH, 'edgecolor', 'none');
view([-90 0]);
lighting gouraud;
material dull;
h1 = light('position',[-100 0 50]); 
h2 = light('position',[-100 0 -50]);
set(gcf,'color','w');
f=gcf;
f.WindowState='maximized';
colormap(cmap_onset);
 plotname = fullfile(savedir, sprintf('%s',model_name), sprintf('onset%s_brain_LateralRight', model_name));
    exportgraphics(f,sprintf('%s.png',plotname),'Resolution',300)

figure('Name', 'Significant Effects Left Medial')
ft_plot_mesh(mesh_left_hemi, 'vertexcolor', vertex_colors_LH,'edgecolor', 'none');
view([90 0]);
lighting gouraud;
material dull;
h1= light('position', [100 0 50]);
h2= light('position', [100 0 -50]);
set(gcf,'color','w');
f=gcf;
f.WindowState='maximized';
colormap(cmap_onset);
 plotname = fullfile(savedir, sprintf('%s',model_name), sprintf('onset%s_brain_MedialLeft', model_name));
    exportgraphics(f,sprintf('%s.png',plotname),'Resolution',300)

figure ('Name', 'Significant Effects Right Medial')
ft_plot_mesh(mesh_right_hemi, 'vertexcolor', vertex_colors_RH, 'edgecolor', 'none');
view([90 0]);
lighting gouraud;
material dull;
h1= light('position', [100 0 50]);
h2= light('position', [100 0 -50]);
set(gcf, 'Color', 'w');
f=gcf;
f.WindowState='maximized';
plotname = fullfile(savedir, sprintf('%s',model_name), sprintf('onset%s_brain_MedialRight', model_name));
exportgraphics(f,sprintf('%s.png',plotname),'Resolution',300)
end