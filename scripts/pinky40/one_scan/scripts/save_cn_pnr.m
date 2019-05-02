%% data before & after running EASE 
Y_cnmf = neuron.normalize_data(Y); 
Y_cnmf = Y_cnmf(:, neuron.frame_range(1):neuron.frame_range(2)); 
Y_res = bsxfun(@minus, Y_cnmf-neuron.b*neuron.f - neuron.A*neuron.C, neuron.b0);

%% summary images before running EASE 
cn_before = cell(1, 3);
pnr_before = cell(1, 3);
Y_cnmf = neuron.reshape(Y_cnmf, 3); 
for m=1:3
    cn_before{m} = correlation_image(squeeze(Y_cnmf(:, :, m, :)));
    pnr_before{m} = pnr_image(squeeze(Y_cnmf(:, :, m, :)));
end

%% summary images after running EASE 
cn_after = cell(1, 3);
pnr_after = cell(1, 3);
Y_res = neuron.reshape(Y_res, 3);
for m=1:3
    cn_after{m} = correlation_image(squeeze(Y_res(:, :, m, :)));
    pnr_after{m} = pnr_image(squeeze(Y_res(:, :, m, :)));
end

%% find EM boundary
ind = neuron.reshape(neuron.spatial_range, 3);
d3 = neuron.options.d3; 
em_ranges = cell(1, d3);
for m=1:3
    [y, x] = find(ind(:, :, m));
    k = convhull(x, y);
    em_ranges{m} = [x(k), y(k)];
end

temp = cell2mat(em_ranges');
xrange = [floor(min(temp(:,1)))-2, ceil(max(temp(:,1)))+2];
yrange = [floor(min(temp(:,2)))-2, ceil(max(temp(:,2)))+2];

%% correlation image before
clear add_scale_bar;
img = cell2mat(cn_before); %#ok<*NASGU>
vlim = [0, 1];
fig_x0 = 10;
fig_y0 = 10;
script_show_image;
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('cn_before_%s.pdf', orientation')));
    export_fig(gcf, fullfile(output_figs, sprintf('cn_before_%s.fig', orientation')));
else
    pause;
end

% correlation image after
img = cell2mat(cn_after);
script_show_image;
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('cn_after_%s.pdf', orientation')));
    export_fig(gcf, fullfile(output_figs, sprintf('cn_after_%s.fig', orientation')));
else
    pause;
end

% generate a colorbar
figure('papersize', [2, 4]);
init_fig;
axes('position', [-0.4, 0.1, 1, 0.8]);
imagesc(rand(100, 0), vlim);
axis off;
colorbar;
colormap jet;
if save_image
    export_fig(gcf, fullfile(output_figs, 'cn_colorbar.pdf'));
    export_fig(gcf, fullfile(output_figs, 'cn_colorbar.fig'));
end

% pnr image before
img = cell2mat(pnr_before);
vlim = [0, 20];
script_show_image;
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('pnr_before_%s.pdf', orientation')));
    export_fig(gcf, fullfile(output_figs, sprintf('pnr_before_%s.fig', orientation')));
else
    pause;
end

% pnr image after
img = cell2mat(pnr_after);
add_scale_bar = true;
pixel_size = ease.range_2p(1) / ease.dims_video(1);
script_show_image;
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('pnr_after_%s.pdf', orientation')));
    export_fig(gcf, fullfile(output_figs, sprintf('pnr_after_%s.fig', orientation')));
else
    pause;
end

% generate a colorbar
figure('papersize', [2, 4]);
init_fig;
axes('position', [-0.4, 0.1, 1, 0.8]);
imagesc(rand(100, 0), vlim);
axis off;
colorbar;
colormap jet;
if save_image
    export_fig(gcf, fullfile(output_figs, 'pnr_colorbar.pdf'));
    export_fig(gcf, fullfile(output_figs, 'pnr_colorbar.fig'));
end