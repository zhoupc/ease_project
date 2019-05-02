
%% show components extracted by
em_ids1 = cell2mat(neuron_ease.match_status.em_ids);
em_ids2 = cell2mat(neuron_cnmf.match_status.em_ids);
em_ids_both = intersect(em_ids1, em_ids2); 
em_ids_ease = setdiff(em_ids1, em_ids2); 
em_ids_cnmf = setdiff(em_ids2, em_ids1); 
t = (1:size(neuron_ease.C, 2))/neuron.Fs;

%%
img = neuron.reshape(neuron.spatial_range, 3);
[img, rot_info] = neuron.rotate_zoomin(img);
[d1, d2, d3] = size(img);
h = 80;
dh = 20;
w = (d2*d3)/d1 * h;

%% visualize neurons that have been assigned to the same EM
tmp_folder = fullfile(output_figs, 'both_extracted');
if ~exist(tmp_folder, 'dir')
    mkdir(tmp_folder);
end
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m);
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id);
    fig_h = 2*(h+dh)* (2+ length(idx_cnmf));
    figure('papersize', [w, fig_h]/100, 'position', [5, 0, w, fig_h]);
    
    % show EM
    img_pi = neuron.rotate_zoomin(neuron_ease.A_em(:, idx_ease), rot_info);
    h0 = fig_h - h- dh;
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_pi, d1, d2*d3));
    axis equal off tight;
    title('EM footprints');
    
    % show correlation image
    h0 = h0 - h- dh;
    img_corr = neuron.rotate_zoomin(neuron_ease.A_corr(:, idx_ease), rot_info);
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_corr, d1, d2*d3));
    axis equal off tight;
    title('correlation image');
    
    %show EASE image
    h0 = h0 - h- dh;
    img_corr = neuron.rotate_zoomin(neuron_ease.A(:, idx_ease), rot_info);
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_corr, d1, d2*d3));
    axis equal off tight;
    if neuron_ease.labels(idx_ease)<0
        title('extracted neuron by EASE', 'color', 'k');
    else
        title('extracted neuron by EASE', 'color', 'k');
    end
    h0 = h0 - h-dh/2;
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    plot(t, neuron_ease.C_raw(idx_ease,:), 'b');
    hold on;
    plot(t, neuron_ease.C(idx_ease,:), 'r');
    axis  tight;
    set(gca, 'xtick', [], 'ytick', []);
    box on;
    
    % show CNMF neurons
    for n=1:length(idx_cnmf)
        h0 = h0 - h- dh;
        img_corr = neuron.rotate_zoomin(neuron_cnmf.A(:, idx_cnmf(n)), rot_info);
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        imagesc(reshape(img_corr, d1, d2*d3));
        axis equal off tight;
        if neuron_cnmf.labels(idx_cnmf(n))<0
            title(sprintf('extracted neuron %d by CNMF',n), 'color', 'k');
        else
            title(sprintf('extracted neuron %d by CNMF',n), 'color', 'k');
        end
        h0 = h0 - h-dh/2;
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        plot(t, neuron_ease.C_raw(idx_ease,:), 'b');
        hold on;
        plot(t, neuron_ease.C(idx_ease,:), 'r');
        axis  tight;
        set(gca,'ytick', []);
        box on;
        set(gca, 'xtick', [], 'ytick', []);
    end
    set(gca, 'xtick', 0:200:t(end));
    
    if save_image
        export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d.pdf', tmp_id)));
        export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d.fig', tmp_id)));
    end
    close;
end

%% visualize neurons that is only visible by EASE
tmp_folder = fullfile(output_figs, 'ease_extracted_only');
if ~exist(tmp_folder, 'dir')
    mkdir(tmp_folder);
end
for m=1:length(em_ids_ease)
    tmp_id = em_ids_ease(m);
    idx_ease = find(em_ids1==tmp_id);
    fig_h = 4*(h+dh)+dh;
    figure('papersize', [w, fig_h]/100, 'position', [5, 0, w, fig_h]);
    
    % show EM
    img_pi = neuron.rotate_zoomin(neuron_ease.A_em(:, idx_ease), rot_info);
    h0 = fig_h - h- dh;
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_pi, d1, d2*d3));
    axis equal off tight;
    title('EM footprints');
    
    % show correlation image
    h0 = h0 - h- dh;
    img_corr = neuron.rotate_zoomin(neuron_ease.A_corr(:, idx_ease), rot_info);
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_corr, d1, d2*d3));
    axis equal off tight;
    title('correlation image');
    
    %show EASE image
    h0 = h0 - h- dh;
    img_corr = neuron.rotate_zoomin(neuron_ease.A(:, idx_ease), rot_info);
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    imagesc(reshape(img_corr, d1, d2*d3));
    axis equal off tight;
    if neuron_ease.labels(idx_ease)<0
        title('extracted neuron by EASE', 'color', 'k');
    else
        title('extracted neuron by EASE', 'color', 'k');
    end
    h0 = h0 - h-dh/2;
    axes('unit', 'pixel', 'position', [5, h0, w, h]);
    plot(t, neuron_ease.C_raw(idx_ease,:), 'b');
    hold on;
    plot(t, neuron_ease.C(idx_ease,:), 'r');
    axis  tight;
    set(gca,'ytick', []);
    box on;
    set(gca, 'xtick', 0:200:t(end));
    xlabel('Time (sec.)');
    
    if save_image
        export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d.pdf', tmp_id)));
        export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d.fig', tmp_id)));
    end
    close;
end

%% visualize neurons that is only visible by CNMF
tmp_folder = fullfile(output_figs, 'cnmf_extracted_only');
if ~exist(tmp_folder, 'dir')
    mkdir(tmp_folder);
end
fig_h = 4*(h+dh)+dh;
figure('papersize', [w, fig_h]/100, 'position', [5, 0, w, fig_h]);

for m=1:length(em_ids_cnmf)
    tmp_id = em_ids_cnmf(m);
    idx_cnmf = find(em_ids2==tmp_id);
    
    for n=1:length(idx_cnmf)
        % show EM
        img_pi = neuron.rotate_zoomin(neuron_cnmf.A_em(:, idx_cnmf(n)), rot_info);
        h0 = fig_h - h- dh;
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        imagesc(reshape(img_pi, d1, d2*d3));
        axis equal off tight;
        title('EM footprints');
        
        % show correlation image
        h0 = h0 - h- dh;
        img_corr = neuron.rotate_zoomin(neuron_cnmf.A_corr(:, idx_cnmf(n)), rot_info);
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        imagesc(reshape(img_corr, d1, d2*d3));
        axis equal off tight;
        title('correlation image');
        
        %show EASE image
        h0 = h0 - h- dh;
        img_corr = neuron.rotate_zoomin(neuron_cnmf.A(:, idx_cnmf(n)), rot_info);
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        imagesc(reshape(img_corr, d1, d2*d3));
        axis equal off tight;
        if neuron_ease.labels(idx_ease)<0
            title('extracted neuron by CNMF', 'color', 'k');
        else
            title('extracted neuron by CNMF', 'color', 'k');
        end
        h0 = h0 - h-dh/2;
        axes('unit', 'pixel', 'position', [5, h0, w, h]);
        plot(t, neuron_cnmf.C_raw(idx_cnmf(n),:), 'b');
        hold on;
        plot(t, neuron_cnmf.C(idx_cnmf(n),:), 'r');
        axis  tight;
        set(gca,'ytick', []);
        box on;
        set(gca, 'xtick', 0:200:t(end));
        xlabel('Time (sec.)');
        
        if save_image
            export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d_%d.pdf', tmp_id, n)));
            export_fig(gcf, fullfile(tmp_folder, sprintf('em_%d_%d.fig', tmp_id, n)));
        end
        clf;
    end
end

%% compare confidence scores 
confidences_ease = neuron_ease.match_status.confidence; 
confidences_cnmf = neuron_cnmf.match_status.confidence; 
figure; hold on; 
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m); 
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id); 
    
    plot(confidences_cnmf(idx_cnmf),confidences_ease(idx_ease),  'ob'); 
end 
plot([0, 6], [0, 6], 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
axis tight; 
title('match confidences'); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'confidences_compare.pdf'));
    export_fig(gcf, fullfile(output_figs, 'confidences_compare.fig'));
end

%% compare correlations 
figure; hold on; 
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m); 
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id); 
    
    pj = neuron_ease.A_em(:, idx_ease); 
    corr_ease = corr(neuron_ease.A(:, idx_ease), pj); 
    corr_cnmf = corr(neuron_cnmf.A(:, idx_cnmf), pj); 
    plot(corr_cnmf, corr_ease,  'ob'); 
end 
plot([0, 1], [0, 1], 'r'); 
xlabel('CNMF'); 
ylabel('EASE'); 
axis tight; 
title('correlation(ai, pi)'); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'correlations_compare.pdf'));
    export_fig(gcf, fullfile(output_figs, 'correlations_compare.fig'));
end

%% compare PNR 
figure; hold on; 
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m); 
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id); 
    
    plot(pnr_cnmf(idx_cnmf), pnr_ease(idx_ease),  'ob'); 
end 
axis tight; 
plot([0, 100], [0, 100], 'r'); 
xlabel('CNMF'); 
ylabel('EASE');
title('PNR'); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'pnr_compare.pdf'));
    export_fig(gcf, fullfile(output_figs, 'pnr_compare.fig'));
end

%% compare snr
figure; hold on; 
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m); 
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id); 
    
    plot(snr_cnmf(idx_cnmf), snr_ease(idx_ease),  'ob'); 
end 
axis tight; 

plot([0, 8], [0, 8], 'r'); 
xlabel('CNMF'); 
ylabel('EASE');
title('SNR'); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'snr_compare.pdf'));
    export_fig(gcf, fullfile(output_figs, 'snr_compare.fig'));
end

%% compare RSS 
figure; hold on; 
for m=1:length(em_ids_both)
    tmp_id = em_ids_both(m); 
    idx_ease = find(em_ids1==tmp_id);
    idx_cnmf = find(em_ids2==tmp_id); 
    
    plot(1-tc_cnmf.rss(idx_cnmf), 1-tc_ease.rss(idx_ease),  'ob'); 
end 
axis tight; 

plot([0, 1], [0, 1], 'r'); 
xlabel('CNMF'); 
ylabel('EASE');
title('SNR'); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'r2_compare.pdf'));
    export_fig(gcf, fullfile(output_figs, 'r2_compare.fig'));
end

%% show distributions 
flag_ease = false(K_ease,1); 
flag_cnmf = false(K_cnmf,1); 
for m=1:length(em_ids_both)
    flag_ease(em_ids1==em_ids_both(m)) = true; 
    flag_cnmf(em_ids2==em_ids_both(m)) = true; 
end 
figure; 
subplot(211); 
hist(tc_ease.rss(~flag_ease), 20); 

subplot(212); 
hist(tc_cnmf.rss(~flag_cnmf), 20); 