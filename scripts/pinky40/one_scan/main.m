%% prepare for running EASE 
if ~exist('ease', 'var')
    fi.usepkg('ease'); 
    run_ease; 
end
addpath('./scripts'); 

ease.scan_id = 1;
ease.block_id = 1; 
ease.use_denoise = false;
ease.create_new = false; 
run_analysis_again = false; 

% prepare for exporting the results 
save_image = true; %#ok<NASGU>
output_figs = fullfile(ease.fig_folder, 'scan_1_raw');
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end

%% run EASE automatically 
if run_analysis_again 
    neuron = ease.get_MF3D(true);
    neuron.frame_range = [201, 9100]; 
else
    neuron = ease.get_MF3D(); 
end

if isempty(neuron.A)
    % configurations of  running EASE 
    ease.options_init.clear_results = true; 
    ease.options_init.save_fig = false; 
    ease.options_init.show_fig = false; 
    ease.options_init.order_statistics = 'l3norm'; 
    ease.options_init.min_similarity = 0.55; 
    neuron.frame_range = [201, 9100];
    niters = 3;
    K_new = [100, 50, 50];
    K_candidate = [2000, 4000, 40000]; 
    nb = [1, 2, 3];  % number of background components
    min_ranks = nan * ones(niters, 1);  % minimum rank of the EM match
    configs = ease.create_running_configs(K_new, K_candidate, nb, min_ranks);
    
    % start running EASE automatically
    tic; 
    neuron.at_ease(configs);
    fprintf('running time: %.2f\n', toc); 
    % running time: 318.61
    ease.save_MF3D(); 
end

%% add some manual intervention 
if run_analysis_again 
    run_interactively;
end


%% summarize results 
fprintf('EASE extracted %d components in total\n', size(neuron.A,2)); 
labels = neuron.labels; 
n_verified = sum(labels<0); 
n_soma = sum(labels==-1); 
n_dendrite = sum(labels==-2); 
fprintf('We manually verified %d matched pairs. \nAmong these matches, there are %d somas and %d dendrites.\n', ...
    n_verified, n_soma, n_dendrite); 

%% save correlation image and pnr image
if ~exist('cn_pnr_saved', 'var') || ~cn_pnr_saved
    save_cn_pnr;
    cn_pnr_saved = true; 
end

%% save example neurons 
if ~exist('examples_saved', 'var') || ~examples_saved
    save_example_neurons;
    cn_pnr_saved = true; 
end
%% cluster neurons according to spatial correlations 
if ~exist('spatial_cluster_saved', 'var') || ~spatial_cluster_saved
    save_spatial_cluster; 
    spatial_cluster_saved = true; 
end 

%% show neurons that are potentially be one neuron  
if ~exist('temporal_cluster_saved', 'var') || ~temporal_cluster_saved
    save_temporal_cluster; 
    temporal_cluster_saved = true; 
end 

%% overlap neuron by grouping them according to confidence values 
if ~exist('order_confidence_saved', 'var') || ~order_confidence_saved 
    save_footprints_order_by_confidence; 
    order_confidence_saved = true; 
end
% %% check the importance of EM-weighting 
% if ~exist('weighting_results', 'file') 
%     weighting_is_important; 
% end 
% 
% %% check the importance of EM info 
% if ~exist('EMinfo_results', 'file') 
%     EMinfo_is_important; 
% end 

%% save a video 
if ~demixing_video_saved
    min_max = [0, 12];
    avi_nm = fullfile(ease.video_folder, ...
        sprintf('demixing_scan_%d_%d_%d.avi', ease.scan_id, ...
        neuron.frame_range(1), neuron.frame_range(2)));
    t_pause = 0;
    Y = ease.load_Y; 
    neuron.showDemixing(Y, min_max, [], avi_nm, t_pause); % ind_neuron, rot_info, fig_visible)
    demixing_video_saved = true; 
end

%% save all neurons
if ~exist('all_neurons_saved', 'var') || ~all_neurons_saved
    save_all_neurons;
    all_neurons_saved = true;
end


