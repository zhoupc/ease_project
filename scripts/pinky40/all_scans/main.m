%% configure parameters
scan_list = 1:4;
scan_to_show = length(scan_list);
ease.block_id = 0;
ease.use_denoise = 0;
run_analysis_again = false;     %

% prepare for exporting the results
save_image = true; %#ok<NASGU>
output_figs = fullfile(ease.fig_folder, sprintf('all_scans_block_%d_raw', ease.block_id));
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end
frame_range = [201, 27300];

% configurations of  running EASE
ease.options_init.clear_results = true;
ease.options_init.save_fig = false;
ease.options_init.show_fig = false;
ease.options_init.min_similarity = 0.55;
niters = 3;
K_new = [100, 50, 50];
K_candidate = [2000, 4000, 40000];
nb = [1, 2, 3];  % number of background components
min_ranks = nan * ones(niters, 1);  % minimum rank of the EM match
configs = ease.create_running_configs(K_new, K_candidate, nb, min_ranks);

%% run analysis
if run_analysis_again
    ease.create_new = true;
    
    %% automatically process all scans
    ease.process_scans(scan_list, configs, frame_range);
    
    %% manual curation
    fprintf('you need to run this step using GUI\n');
    ease.startGUI();
    fprintf('type ENTER to continue: ');
    pause();
    
    %% summarize results and create a white list
    results = ease.combine_results(scan_list);
    count_verified = sum(results.labels<0, 1);
    ind = (count_verified>0);
    white_list = results.em_ids(ind);
    
    %% initialize neurons given the white list
    ease.options_init.clear_results = false;
    ease.options_init.show_fig = true;
    configs_res = ease.create_running_configs(length(white_list),...
        length(segment_ids), 3, nan, [], white_list);
    
    % process all scans
    ease.process_scans(scan_list, configs_res);
    
    % compute tuning curves
    shuffle = [];
    ease.compute_tuning_curves(scan_list, shuffle);
end

%% compute tuning curves
[results, neurons_all] = ease.combine_results(scan_list);

em_id = 23683777; 
example_idx = find(results.em_ids==em_id); 
visualize_one_neuron; 
%% order neurons according to confidence values
confidences = max(results.confidences, [], 1);
[~, idx] = sort(confidences, 'descend');
for m=1:length(idx)
    example_idx = idx(m);
    close all;
    disp([m, example_idx]);
    visualize_one_neuron;
end
