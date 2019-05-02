%% prepare for running EASE 
if ~exist('ease', 'var')
    fi.usepkg('ease'); 
    run_ease; 
end

ease.scan_id = 1;
ease.block_id = 1; 
ease.use_denoise = false;
ease.create_new = false; 

% prepare for exporting the results 
save_image = true;
output_figs = fullfile(ease.fig_folder, 'initialization');
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end

%% load video data into the memory
neuron = ease.get_MF3D(); 
ind_frame = 201:9100; 
neuron.frame_range = ind_frame([1, end]); 
Y = ease.load_Y(); 
[Y_cnmf, Y_sn] = neuron.normalize_data(Y);

[Aem, segment_ids] = ease.load_Aem(); 

%% load test data and run test analysis   
load test_data; 
K_test = size(A_2p, 2); 
example_id = 3; 
col = parula(7); 
markers = 'ovds>'; 
pixel_size = ease.range_2p(1) / ease.dims_video(1); 

% set ground truth 
ai_test = A_2p(:, example_id);
ai_em_test = A_em(:, example_id); 
ci_test = C_test(example_id,:) * 0.6; 
ci_test = ci_test - mean(ci_test); 

% get the selected EM neuron and extract video data 
ai_em = neuron.reshape(ai_em_test, 3) / sum(ai_em_test(:).^2);
% get pixels inside the mask and outside of the mask
ai_mask = imdilate(ai_em>0, strel('square', 1));  
ind_in = neuron.reshape(ai_mask, 3);
ind_out = xor(ind_in, imdilate(ai_mask, strel('square', 8)));
Yin = Y_cnmf(ind_in(:), ind_frame); %-neuron.b(ind_in(:), :)*neuron.f(ind_frame) ;
Yout = Y_cnmf(ind_out(:), ind_frame); % -neuron.b(ind_out(:), :)*neuron.f(ind_frame);
Yin = bsxfun(@minus, Yin, mean(Yin, 2));
Yout = bsxfun(@minus, Yout, mean(Yout, 2));

% results on the original SNR level
Yin_new = Yin + ai_test(ind_in(:))*ci_test; 
ci0 = reshape(ai_em(ind_in(:)), 1, []) * Yin_new ;
rank_out = 100; 

% extract spatial &  temporal activity
ai_new = double(ind_in);
% ease method, wrss 
ai_tf = double(ind_in);     
[ai_tf(ind_in), ci_tf] = initialize_ac_tf(Yin_new, Yout, ai_em(ind_in(:)), ci0, rank_out);

% ease method, rss 
ai_tf_uniform = double(ind_in);  
[ai_tf_uniform(ind_in), ci_tf_uniform] = initialize_ac_tf(Yin_new, Yout, [], ci0, rank_out);

% semi-nmf method, wrss 
ai_iter = double(ind_in);   % use semi-NMF method
[ai_iter(ind_in), ci_iter] = initialize_ac_seminmf(Yin_new, ai_em(ind_in(:)), ci0);

% semi-nmf method, rss 
ai_iter_uniform = double(ind_in);   % use semi-NMF method
[ai_iter_uniform(ind_in), ci_iter_uniform] = initialize_ac_seminmf(Yin_new, [], ci0);

% put all figures into one image  
tmp_plane_id = 1;   % only show one plane  
img = zeros(numel(ai_new), 5); 
img(:, 1) = ai_test(:)/ norm(ai_test,2); %sum(ci_test); 
img(:, 2) = ai_iter_uniform(:) / norm(ai_iter_uniform(:), 2); %*sum(ci_iter_uniform); 
img(:, 3) = ai_iter(:) / norm(ai_iter(:), 2); %* sum(ci_iter); 
img(:, 4) = ai_tf_uniform(:) / norm(ai_tf_uniform(:),2); %* sum(ci_tf_uniform); 
img(:, 5) = ai_tf(:) / norm(ai_tf(:), 2); %* sum(ci_tf);
img = [ai_em_test(:)/norm(ai_em_test(:),2), img]; 
img = neuron.reshape(img, 3); 
[rlim , clim] = determine_bounding_box(sum(img, 4)); 
rrange = rlim(1):rlim(2);
crange = clim(1):clim(2);
img = squeeze(img(rrange, crange, tmp_plane_id, :)); 
neuron.showImage(img, 'horizontal', [], pixel_size); 

if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('results_ai_%d.pdf', example_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('results_ai_%d.fig', example_id)));
end 

% temporal traces 
ci_true = ci_test / std(ci_test); 
T = length(ci_true); 
tt = (1:T) / neuron.Fs; 
trange = [100, 350]; 
ind = (tt>trange(1)) & (tt<trange(2));  % choose bins to display 
tt = tt(ind); 
scale_bar = diff(trange) / 10; 
tmp_amp = ceil(max(ci_true(ind))); 

% temporal 
figure('papersize', [9, 3]); 
init_fig;
axes('position', [0.005, 0.005, 0.994, 0.99]); hold on; 
set(gca, 'fontsize', 13);
ci_tf_denoised = deconvolveCa(ci_tf, neuron.options.deconv_options); %, 'smin', -3);
tmp_C = {ci0, ci_iter_uniform, ci_iter, ci_tf_uniform, ci_tf, ci_tf_denoised};
tmp_str = {'mean', 'semi-NMF,RSS', 'semi-NMF,wRSS', 'EASE,RSS', 'EASE,wRSS', ...
     'EASE,RSS,denoised', 'true'};
for m=1:6
    if m==1
        tmp_col = 'g'; 
            temp = tmp_C{m} / std(ci_test); 
    else
        tmp_col = col(m-1,:);
            temp = normalize(tmp_C{m}, 'zscore'); 
    end
    p = plot(tt, tmp_amp*(m-1) + temp(ind), '-', 'color', tmp_col);
    text(trange(end)-20, tmp_amp*m-4, sprintf('%.3f', ...
         corr(tmp_C{m}(:), ci_test(:))), 'color', tmp_col, ...
        'fontsize', 13); 
end 
for m=1:6
   plot(tt, tmp_amp*(m-1) + ci_true(ind), '-.r', 'linewidth', 1); 
end

[tmp_h, tmp_labels] = legend('mean', 'semi-NMF,RSS', 'semi-NMF,wRSS', 'EASE,RSS', 'EASE,wRSS', ...
    'EASE,RSS,denoised', 'true', 'orientation', 'horizontal'); 
for m=1:7 
    tmp_labels(6+2*m).XData = tmp_labels(6+2*m).XData - [0.02*(m-1), 0.02*m]; 
    tmp_labels(m).Position = tmp_labels(m).Position - [0.02*m, 0, 0]; 
end 
temp = [0.14, 0.9, 0.9, 0.1]; 
set(tmp_h, 'position', temp, 'Box', 'off', 'fontsize', 13); 

set(gca, 'xtick', []); 
set(gca, 'ytick', []); 
axis tight; 
xlim(trange); 
ylim([-tmp_amp/2, tmp_amp*6.5]); 
box on; 
% scale bar 
plot(tt(end) + [-30, -10], [1, 1]*(-tmp_amp/3), 'k', 'linewidth', 5); 

if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('results_ci_%d.pdf', example_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('results_ci_%d.fig', example_id)));
end 

% try different SNR levels 
snr_factors = 0.1:0.1:2; 
num_factors = length(snr_factors); 
corr_ai = zeros(5, num_factors); 
corr_ci = zeros(5, num_factors); 

for k_factor = 1:num_factors
    Yin_new = Yin + ai_test(ind_in(:))*ci_test * snr_factors(k_factor); 
    ci0 = reshape(ai_em(ind_in(:)), 1, []) * Yin_new ;

    %% extract spatial &  temporal activity
    ai_new = double(ind_in);
    ai_tf = double(ind_in);     % use TF method
    ai_tf_uniform = double(ind_in);     % use TF method
    ai_iter = double(ind_in);   % use semi-NMF method
    ai_iter_uniform = double(ind_in);   % use semi-NMF method
    [ai_tf(ind_in), ci_tf] = initialize_ac_tf(Yin_new, Yout, ai_em(ind_in(:)), ci0, 100);
    [ai_tf_uniform(ind_in), ci_tf_uniform] = initialize_ac_tf(Yin_new, Yout, [], ci0, 100);
    [ai_iter(ind_in), ci_iter] = initialize_ac_seminmf(Yin_new, ai_em(ind_in(:)), ci0);
    [ai_iter_uniform(ind_in), ci_iter_uniform] = initialize_ac_seminmf(Yin_new, [], ci0);
    
    %% compute correlation coefficients
    corr_ai(4, k_factor) = corr(ai_tf(:), ai_test(:));
    corr_ai(3, k_factor) = corr(ai_tf_uniform(:), ai_test(:));
    corr_ai(2, k_factor) = corr(ai_iter(:), ai_test(:));
    corr_ai(1, k_factor) = corr(ai_iter_uniform(:), ai_test(:));
    
    ci_tf_denoising = deconvolveCa(ci_tf, neuron.options.deconv_options);
    corr_ci(5, k_factor) = corr(ci_tf_denoising(:), ci_test(:));
    corr_ci(4, k_factor) = corr(ci_tf(:), ci_test(:));
    corr_ci(3, k_factor) = corr(ci_tf_uniform(:), ci_test(:));
    corr_ci(2, k_factor) = corr(ci_iter(:), ci_test(:));
    corr_ci(1, k_factor) = corr(ci_iter_uniform(:), ci_test(:));
    
end

% spatial correlation 
figure('papersize', [7,5]); 
init_fig; 
axes('position', [0.11, 0.13, 0.86, 0.85]); hold on; 
set(gca, 'fontsize', 17); 
for m=1:4
    plot(snr_factors, corr_ai(m,:), 'marker', markers(m),...
        'markersize', 8, 'color', col(m,:)); 
end 
box on; 
axis([0, snr_factors(end)+0.1, 0.0, 1]); 
set(gca, 'xtick', 0:0.5:snr_factors(end)); 
set(gca, 'ytick', 0:0.2:1); 
xlabel('SNR factor'); %, 'interpreter', 'latex');
ylabel('corr(true, initialized)'); %, 'interpreter', 'latex'); 
% ylabel('corr($\hat{a}_i, a_i$)', 'interpreter', 'latex'); 
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('ai_corr_%d.pdf', example_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('ai_corr_%d.fig', example_id)));
end 

% temporal 
figure('papersize', [7,5]); 
init_fig; 
axes('position', [0.11, 0.13, 0.86, 0.85]); hold on; 
set(gca, 'fontsize', 17); 
for m=1:5
    plot(snr_factors, corr_ci(m,:), 'marker', markers(m), ...
        'markersize', 8, 'color', col(m,:)); 
end 
box on; 
axis([0, snr_factors(end)+0.05, -0.02, 1]); 
set(gca, 'xtick', 0:0.5:snr_factors(end)); 
set(gca, 'ytick', 0:0.2:1);
xlabel('SNR factor'); %, 'interpreter', 'latex');
% ylabel('corr(true, inferred)'); %, 'interpreter', 'latex'); 
% xlabel('SNR factor', 'interpreter', 'latex');
% ylabel('corr($\hat{c}_i, c_i$)', 'interpreter', 'latex'); 


tmp_h = legend('semi-NMF,RSS', 'semi-NMF,wRSS', 'EASE,RSS', 'EASE,wRSS', ...
    'EASE,wRSS,denoised', 'orientation', 'vertical', 'fontsize', 14); 
temp = get(tmp_h, 'position');
temp(2) = 0.45; 
temp(1) = 0.93 - temp(3); 
set(tmp_h, 'position', temp, 'fontsize', 14); 
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('ci_corr_%d.pdf', example_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('ci_corr_%d.fig', example_id)));
end 
