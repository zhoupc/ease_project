%% configure parameters
ease.scan_id = 1; 
ease.block_id = 0;
ease.use_denoise =false;
save_image = true;
run_analysis_again = false;
ease.create_new = false;
output_figs = fullfile(ease.fig_folder, sprintf('compare_cnmf_scan_%d_new', ease.scan_id));
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end

% load stimuli conditions 
ease.load_stimuli(); 
plane_id = 1:3; 

% load data 
cb_btn_load_ca;

%% load EASE results
if ~exist('ease_results.mat', 'file')
    % fetch results 
    neuron_ease = ease.get_MF3D();
    
    % update teporal traces again; this time, no weights from EM 
    neuron_ease.update_temporal(Y, false, false, true);
    
    % evaluate matching confidences 
    neuron_ease.evaluate_matching_confidence(Y, Aem, segment_ids); 
    
    % verify matches 
    ease.startGUI();
    fprintf('do manual verification\n');
    neuron = neuron_ease.copy();
    pause; 
    neuron_ease = neuron.copy(); 

    % compute tuning curves 
    neuron_ease.compute_tuning_curve(ease.conditions(ease.scan_id, 201:end));
    neuron_ease.boostrap_tuning_curve(nshuffle);

    % copy results 
    K_ease = size(neuron_ease.A, 2);
    A_ease = neuron_ease.reshape(neuron_ease.A, 3);
    A_ease = reshape(A_ease(:, :, plane_id, :), [], K_ease);
    C_ease = neuron_ease.C;
    Craw_ease = neuron_ease.C_raw;     
    tc_ease = neuron_ease.tuning_curve;
    scores_ease = neuron_ease.match_status.scores;
    confidences_ease =neuron_ease.match_status.confidence; 
    em_ids_ease = cell2mat(neuron_ease.match_status.em_ids); 
    labels_ease = neuron_ease.scores; 
    save ease_results.mat neuron_ease  K_ease A_ease C_ease Craw_ease ...
        tc_ease scores_ease;
else
    load ease_results.mat;
end
%% load CNMF results
nshuffle = []; 
if ~exist('cnmf_results.mat', 'file')
    % create a wrapper for storing CNMF result
    neuron_cnmf = neuron.copy();
    
    % fetch results from Ding 
    rlt = matfile('rlt_update.mat');
    temp = neuron.reshape(neuron.spatial_range, 3);
    tmp_ind = neuron.reshape(neuron.spatial_range, 3);
    tmp_ind(13:52, 1:128, plane_id) = true;
    d = numel(tmp_ind); 
    K_cnmf = size(rlt.a, 2);
    neuron_cnmf.A = zeros(d, K_cnmf);
    neuron_cnmf.A(tmp_ind, :) = rlt.a;
    neuron_cnmf.C = rlt.c';
    neuron_cnmf.C_raw = rlt.c';
    neuron_cnmf.b = zeros(d, size(rlt.fb,2));
    neuron_cnmf.b(tmp_ind,:) = rlt.fb;
    neuron_cnmf.options.nb = size(neuron_cnmf.b, 2);
    neuron_cnmf.f = rlt.ff';
    neuron_cnmf.b0 = zeros(d,1);
    neuron_cnmf.b0(tmp_ind, 1) = rlt.b;
    
    % update teporal traces again
    neuron_cnmf.update_temporal(Y, false, false, true);
    
    % evaluate matching confidences 
    neuron_cnmf.evaluate_matching_confidence(Y, Aem, segment_ids); 
    
    % match to the EM components according to their matching scores 
    [vmax, idx]= max(neuron_cnmf.scores, [],2); 
    em_ids = segment_ids(idx);
    neuron_cnmf.A_em = Aem(:, idx); 
    neuron_cnmf.match_status.em_ids = num2cell(em_ids); 
    neuron_cnmf.match_status.status = zeros(size(em_ids)); 
    neuron_cnmf.match_status.scores = full(vmax); 
    v = sort(neuron_cnmf.scores, 2); 
    neuron_cnmf.match_status.confidence = full(v(:, end)./v(:, end-1));  
    neuron_cnmf.labels = zeros(size(em_ids)); 
    
    % do manual verification
    ease.startGUI;
    neuron = neuron_cnmf.copy();
    fprintf('manually rematch components.\n');
    pause;
    neuron_cnmf = neuron.copy(); 
    
    % compute tuning curves 
    neuron_cnmf.compute_tuning_curve(ease.conditions(ease.scan_id, 201:end));
    neuron_cnmf.boostrap_tuning_curve(nshuffle);

    % copy results 
    K_cnmf = size(neuron_cnmf.A, 2);
    A_cnmf = neuron_cnmf.reshape(neuron_cnmf.A, 3);
    A_cnmf = reshape(A_cnmf(:, :, plane_id, :), [], K_cnmf);
    C_cnmf = neuron_cnmf.C;
    Craw_cnmf = neuron_cnmf.C_raw;     
    tc_cnmf = neuron_cnmf.tuning_curve;
    scores_cnmf = max(neuron_cnmf.scores, [], 1);
    labels_cnmf = neuron_cnmf.labels; 
    
    % save the results
    save cnmf_results.mat neuron_cnmf  K_cnmf A_cnmf C_cnmf Craw_cnmf ...
        tc_cnmf scores_cnmf;
else
    load cnmf_results.mat;
end

%%
pnr_cnmf = max(Craw_cnmf, [], 2) ./ GetSn(Craw_cnmf);
pnr_ease = max(Craw_ease, [], 2) ./ GetSn(Craw_ease);

snr_cnmf = std(C_cnmf, 0, 2) ./ std(Craw_cnmf-C_cnmf, 0, 2);
snr_ease = std(C_ease, 0, 2) ./ std(Craw_ease-C_ease, 0, 2);

l3norm_cnmf = nthroot(sum(A_cnmf.^3, 1)' .* sum(bsxfun(@minus, Craw_cnmf, mean(Craw_cnmf, 2)).^3, 2), 3); 
l3norm_ease = nthroot(sum(A_ease.^3, 1)' .* sum(bsxfun(@minus, Craw_ease, mean(Craw_ease, 2)).^3, 2), 3); 

spatial_range = neuron_cnmf.spatial_range; 
A_corr = corr(A_cnmf(spatial_range, :), A_ease(spatial_range, :)); 
C_corr = corr(C_cnmf', C_ease');

C_corr(isnan(A_corr)) = -inf; 
A_corr(isnan(A_corr)) = -inf; 

% order 
K_cnmf = size(A_cnmf, 2); 
K_ease = size(A_ease, 2); 

[vmax, ~] = max(C_corr, [], 1);
[~, idx_ease] = sort(vmax, 'descend');
idx_cnmf = zeros(K_cnmf, 1); 
for m=1:K_cnmf
    temp = C_corr(:, idx_ease(m)); 
    temp(idx_cnmf(idx_cnmf>0)) = -inf; 
    [~, idx_cnmf(m)] = max(temp); 
end 

%% 
figure('papersize', [4.4, 2]);
init_fig; 
axes('position', [0.02, 0.13, 0.9, 0.86]); 
imagesc(C_corr(idx_cnmf, idx_ease)); 
tmp_pos = get(gca, 'position'); 
colorbar; 
set(gca, 'position', tmp_pos); 
axis equal  tight; 
set(gca, 'xtick', 0:20:180); 
set(gca, 'xticklabel', []); 
set(gca, 'yticklabel', []); 
ylabel('CNMF'); 
xlabel('EASE'); 
% title('temporal correlation'); 
set(gca, 'fontsize', 13); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'temporal_corr.pdf'));
    export_fig(gcf, fullfile(output_figs, 'temporal_corr.fig'));
end
%%
figure('papersize', [4.4, 2]);
init_fig; 
axes('position', [0.02, 0.13, 0.9, 0.86]); 
imagesc(A_corr(idx_cnmf, idx_ease)); 
tmp_pos = get(gca, 'position'); 
colorbar; 
set(gca, 'position', tmp_pos); 
axis equal  tight; 
set(gca, 'xtick', 0:20:180); 
set(gca, 'xticklabel', []); 
set(gca, 'yticklabel', []); 
% ylabel('CNMF'); 
xlabel('EASE'); 
% title('spatial correlation'); 
set(gca, 'fontsize', 13); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial_corr.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_corr.fig'));
end

%% sorted PNR 
idx_sort_ease = idx_ease; 
idx_sort_cnmf = idx_cnmf; 
% [~, idx_sort_ease] = sort(pnr_ease, 'descend'); 
% [~, idx_sort_cnmf] = sort(pnr_cnmf, 'descend'); 

figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
plot(pnr_cnmf(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1);
hold on; 
plot(pnr_ease(idx_sort_ease), 'sr', 'markersize', 5, 'linewidth', 1);

% plot(pnr_ease(idx_ease), 'r'); 
axis tight; box on; 
legend('CNMF', 'EASE'); 
% ylabel('PNR'); 
xlim([0, K_ease]); 
ylim([0,1].*get(gca, 'ylim')); 
set(gca, 'xticklabel', []); 
% xlabel('component #'); 

set(gca, 'fontsize', 14); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_pnr.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_pnr.fig'));
end

%% sorted SNR
idx_sort_ease = idx_ease; 
idx_sort_cnmf = idx_cnmf; 
% [~, idx_sort_ease] = sort(snr_ease, 'descend'); 
% [~, idx_sort_cnmf] = sort(snr_cnmf, 'descend'); 

figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
plot(snr_cnmf(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1);
hold on; 
plot(snr_ease(idx_sort_ease), 'sr', 'markersize', 5, 'linewidth', 1);

% plot(snr_ease(idx_ease), 'r'); 
axis tight; box on; 
legend('CNMF', 'EASE'); 
% ylabel('SNR'); 
xlim([0, K_ease]); 
ylim([0,1].*get(gca, 'ylim')); 

set(gca, 'xticklabel', []); 
% xlabel('component #'); 

set(gca, 'fontsize', 14); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_snr.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_snr.fig'));
end
%% sorted l3 norm 
idx_sort_ease = idx_ease; 
idx_sort_cnmf = idx_cnmf; 
% [~, idx_sort_ease] = sort(l3norm_ease, 'descend'); 
% [~, idx_sort_cnmf] = sort(l3norm_cnmf, 'descend'); 

figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 

plot(l3norm_cnmf(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1);
hold on; 
plot(l3norm_ease(idx_sort_ease), 'sr', 'markersize', 5, 'linewidth', 1);

axis tight; box on; 
legend('CNMF', 'EASE'); 
xlim([0, K_ease]); 

set(gca, 'xticklabel', []); 

set(gca, 'fontsize', 14); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_l3norm.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_l3norm.fig'));
end

%% sorted RSS 
figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
% plot(1-sort(tc_cnmf.rss, 'ascend'), 'o', 'markersize', 5, 'linewidth', 1); 
% hold on; 
% plot(1-sort(tc_ease.rss, 'ascend'), 'sr', 'markersize', 5, 'linewidth', 1); 
plot(1-tc_cnmf.rss(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1); 
hold on; 
plot(1-tc_ease.rss(idx_sort_ease),  'sr', 'markersize', 5, 'linewidth', 1); 
axis tight; box on;
xlim([0, K_ease]); 
% legend('CNMF', 'EASE'); 
% ylabel('RSS'); 
xlabel('rank'); 

set(gca, 'fontsize', 14); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_rss.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_rss.fig'));
end

%%
%% sorted pvals
figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
% plot(1-sort(tc_cnmf.rss, 'ascend'), 'o', 'markersize', 5, 'linewidth', 1); 
% hold on; 
% plot(1-sort(tc_ease.rss, 'ascend'), 'sr', 'markersize', 5, 'linewidth', 1); 
plot(tc_cnmf.pvals(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1); 
hold on; 
plot(tc_ease.pvals(idx_sort_ease),  'sr', 'markersize', 5, 'linewidth', 1); 
axis tight; box on;
xlim([0, K_ease]); 
% legend('CNMF', 'EASE'); 
% ylabel('RSS'); 
xlabel('rank'); 
plot([0, K_ease], [0.05, 0.05], 'color', [1,1,1]*0.6, 'linewidth', 1); 
set(gca, 'fontsize', 14); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_pvals.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_pvals.fig'));
end
%% sorted matching score
figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
% plot(pnr_cnmf(idx_cnmf), ''); 
plot(scores_cnmf(idx_cnmf), 'o', 'markersize', 5, 'linewidth', 1); 
hold on; 
plot(scores_ease(idx_ease), 'sr', 'markersize', 5, 'linewidth', 1); 
% plot(pnr_ease(idx_ease), 'r'); 
axis tight; box on; 
xlim([0, K_ease]); 
legend('CNMF', 'EASE'); 
set(gca, 'xticklabel', []); 

set(gca, 'fontsize', 14); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_score.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_score.fig'));
end

figure('papersize', [6.8, 2]);
init_fig; 
axes('position', [0.1, 0.28, 0.85, 0.7]); 
% plot(pnr_cnmf(idx_cnmf), ''); 
plot(scores_cnmf(idx_sort_cnmf), 'o', 'markersize', 5, 'linewidth', 1); 
hold on; 
plot(scores_ease(idx_sort_ease), 'sr', 'markersize', 5, 'linewidth', 1); 
% plot(pnr_ease(idx_ease), 'r'); 
axis tight; box on; 
xlim([0, K_ease]); 
% legend('CNMF', 'EASE'); 
% ylabel('match score'); 
set(gca, 'xticklabel', []); 
% xlabel('component #'); 

set(gca, 'fontsize', 14); 

if save_image
    export_fig(gcf, fullfile(output_figs, 'sorted_score.pdf'));
    export_fig(gcf, fullfile(output_figs, 'sorted_score.fig'));
end
%%
row_max = max(C_corr, [], 2);
col_max = max(C_corr, [], 1); 
ind = bsxfun(@eq, C_corr, row_max) & bsxfun(@eq, C_corr, col_max) & (C_corr>0); 
[r, c] = find(ind); 
figure; 
plot(pnr_cnmf(r),pnr_ease(c),  'o'); 
hold on; 
plot( pnr_ease(c),pnr_ease(c), 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
title('PNR'); 

%%
row_max = max(C_corr, [], 2);
col_max = max(C_corr, [], 1); 
ind = bsxfun(@eq, C_corr, row_max) & bsxfun(@eq, C_corr, col_max) & (C_corr>0); 
[r, c] = find(ind); 
figure; 
plot(snr_cnmf(r),snr_ease(c),  'o'); 
hold on; 
plot(snr_ease(c),snr_ease(c), 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
title('SNR'); 

%% 
% ind = (A_corr>0.7) & (C_corr>0.8); 
% [r, c] = find(ind); 
r = idx_cnmf(1:10); 
c = idx_ease(1:10); 
figure; 
plot(snr_cnmf(r),snr_ease(c),  'o'); 
hold on; 
plot(snr_ease(c),snr_ease(c), 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
title('SNR'); 

%%
r = idx_cnmf(1:10); 
c = idx_ease(1:10); 
figure; 
plot(pnr_cnmf(r), pnr_ease(c),  'o'); 
hold on; 
plot(pnr_ease(c), pnr_ease(c), 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
title('SNR'); 