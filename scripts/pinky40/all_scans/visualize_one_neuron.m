%% visualize spatial footprints  
FOV = ease.FOV; 
neuron = neurons_all{1}; 
d1 = neuron.options.d1; 
d2 = neuron.options.d2; 
d3 = neuron.options.d3; 
spatial_res = ease.range_2p(1:2) ./ ease.dims_video; 
z_zoomin = 2;
view_angle = [30, 28]; 
zvals = reshape(ease.video_zvals_updated(scan_list,:)', 1, []); 
col_select = cool(5); 

% load saved results and cound neurons repeated in
V_2p = nan(d1, d2, d3, scan_to_show);
V_em = V_2p;
C_raw = cell(scan_to_show, 1);
C = cell(scan_to_show, 1);
S = cell(scan_to_show, 1);
TC = cell(scan_to_show, 1); % tuning curve

em_id = results.em_ids(example_idx); 

idx_each_scan = zeros(scan_to_show, 1);
C_raw = results.Craw(:, :, example_idx); 
C = results.C(:, :, example_idx); 
S = results.S(:, :, example_idx); 
TC = results.y(:, :, example_idx); 
V_2p = reshape(full(results.A(:, example_idx)), d1, d2, []); 
V_em = reshape(full(results.A(:, example_idx)), d1, d2, []); 
V_2p(V_2p==0) = nan; 
V_em(V_em==0) = nan; 
V_2p = reshape(V_2p, d1, d2, []);
V_em = reshape(V_em, d1, d2, []);

fig_size = [5.5, 7]; 

%% show EM meshes
key = struct('segmentation', ease.em_segmentation, 'segment_id', em_id);
% fetch EM meshes
[vertices, faces] = fetchn(ease.rel_mesh & key, 'vertices', 'triangles');
figure('papersize', fig_size, 'name', 'EM meshes');
fig_x0 = 10;
fig_y0 = 400;
init_fig;
hold on; 
for m=1:length(faces)
    vert = bsxfun(@plus, vertices{m} * ease.em_scale_factor * ...
        ease.transformation.A_convert, ease.transformation.offset);
    k_vert = size(vert, 1);
    trisurf(faces{m}+1, vert(:,1)-(FOV(3)-3)*spatial_res(2),...
        ease.range_2p(2)-vert(:,2)-(FOV(1)-3)*spatial_res(2), ...
        (ease.range_2p(3)-vert(:,3))*z_zoomin, ...
        'edgecolor', 'none', 'facecolor', [1, 0.7, 0]);
end
hold on; 

axis equal tight;

% determine xlim and ylim 
xrange =  get(gca, 'xlim')+[-2, 2]; 
yrange =  get(gca, 'ylim')+[-2,2]; 
zrange = get(gca, 'zlim')+[-2,2]; 

[xx, yy] = meshgrid((1:d2)*spatial_res(2), (1:d1)*spatial_res(1)); 
for m=1:(scan_to_show*d3)
    surf(xx, yy, ones(size(xx))*zvals(m)*z_zoomin, nan(size(xx)),...
        'edgecolor', [1,1,1]*0.9, 'edgealpha', 0.5);
    alpha(1);
    hold on;
    
    if mod(m, d3)==0
        plot3(ones(3,1)*xrange(1), ...
            ones(3,1)*yrange(1), zvals(m+(-2:0))*z_zoomin, 'color', ...
            col_select(m/d3, :), 'linewidth', 5);
    end
end
set(gca, 'xtick', (floor(xrange(1)/20):ceil(xrange(2)/20))*20);
set(gca, 'ytick', (floor(yrange(1)/20):ceil(yrange(2)/20))*20);
set(gca, 'ztick', (60:20:180)*z_zoomin);
%     set(gca, 'zticklabel', '');
view(view_angle);
caxis([0, max(V_2p(:))]);

temp = colormap('hot');
colormap(flipud(temp));
grid off;
box on;
set(gca, 'position', [0.16, 0.1, 0.78, 0.88], 'fontsize', 14);
set(gca, 'xlim', xrange, 'ylim', yrange, ...
    'zlim', zrange, ...
    'boxstyle', 'full');
set(gca, 'zticklabel', get(gca, 'ztick')/z_zoomin);

xlabel('X (um)', 'rotation', -20);
ylabel('Y (um)', 'rotation', 30);
zlabel('Z (um)');
ax_pos = get(gca, 'position'); 

if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_meshes_3d.fig', em_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_meshes_3d.png', em_id)));
end

%% show EM footprints 
figure('papersize', fig_size, 'name', '2P footprints');
fig_x0 = 510;
fig_y0 = 400;
init_fig;
for m=1:(scan_to_show*d3)
    surf(xx, yy, ones(size(xx))*zvals(m)*z_zoomin, 4*(V_em(:, :, m)),...
        'edgecolor', [1,1,1]*0.9, 'edgealpha', 0.5);
    alpha(1);
    hold on;
end
axis equal tight;
set(gca, 'xtick', (floor(xrange(1)/20):ceil(xrange(2)/20))*20); 
set(gca, 'ytick', (floor(yrange(1)/20):ceil(yrange(2)/20))*20); 
set(gca, 'ztick', (60:20:180)*z_zoomin);
set(gca, 'zticklabel', '');
view(view_angle);
caxis([0, max(V_em(:))]);

temp = colormap('hot');
colormap(flipud(temp)); 
grid off;
box on;
set(gca, 'position', [0.16, 0.1, 0.78, 0.88], 'fontsize', 14);
set(gca, 'xlim', xrange, 'ylim', yrange, ...
    'zlim', zrange, ...
    'boxstyle', 'full');
xlabel('X (um)', 'rotation', -10);
ylabel('Y (um)', 'rotation', 45);
set(gca, 'position', ax_pos); 

if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_em_3d.fig', em_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_em_3d.png', em_id)));
end

%% show 2P footprints
figure('papersize', fig_size, 'name', '2P footprints');
fig_x0 = 1010;
fig_y0 = 400;
init_fig;

for m=1:(scan_to_show*d3)
    surf(xx, yy, ones(size(xx))*zvals(m)*z_zoomin, 4*(V_2p(:, :, m)),...
        'edgecolor', [1,1,1]*0.9, 'edgealpha', 0.5);
    alpha(1);
    hold on;
end

axis equal tight;
set(gca, 'xtick', (floor(xrange(1)/20):ceil(xrange(2)/20))*20); 
set(gca, 'ytick', (floor(yrange(1)/20):ceil(yrange(2)/20))*20); 
set(gca, 'ztick', (60:20:180)*z_zoomin);
set(gca, 'zticklabel', '');
view(view_angle);
caxis([0, max(V_em(:))]);

temp = colormap('hot');
colormap(flipud(temp)); 
grid off;
box on;
set(gca, 'position', [0.16, 0.1, 0.78, 0.88], 'fontsize', 14);
set(gca, 'xlim', xrange, 'ylim', yrange, ...
    'zlim', zrange, ...
    'boxstyle', 'full');
xlabel('X (um)', 'rotation', -10);
ylabel('Y (um)', 'rotation', 45);
set(gca, 'position', ax_pos); 
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_2p_3d.fig', em_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_2p_3d.png', em_id)));
end

%% show temporal traces
figure('papersize', [8, 4], 'name', 'Temporal traces');
fig_x0 = 100;
fig_y0 = 0;
init_fig;
clear axes;
axes('position', [0.05, 0.01, 0.94, 0.94]);
hold on;
T = size(neuron.C, 2);
tt = (1:T) / neuron.Fs;
trange = [400, 500];
for m=1:scan_to_show
    y_raw = C_raw(:, m); 
    y = C(:, m);
    v_norm = 15; 
    if ~isempty(y_raw)
        plot(tt, y_raw/v_norm+(scan_to_show-m), '-.', 'linewidth', 1, 'color', [1,1,1]*0.3);
        plot(tt, y/v_norm+(scan_to_show-m), 'linewidth', 2, 'color', col_select(m,:));
    end
end
axis tight; 
set(gca, 'xtick', [], 'ytick', 0:4, 'yticklabel', 1:5);
set(gca, 'fontsize', 12);
% scale bar
plot(trange(end) + [-15, -5], [1, 1]*(-0.3), 'k', 'linewidth', 4);
axis tight; 
xlim(trange);
ylim(get(gca, 'ylim')+[-0.1, 0.2]);
box on;
set(gca, 'ytick', (1:scan_to_show)-1);
set(gca, 'fontsize', 17); 
if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_traces', em_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_traces.pdf', em_id)));
end

%% show tuning curves 
figure('papersize', [3, 4], 'name', 'Temporal traces');
fig_x0 = 100;
fig_y0 = 0;
init_fig;
clear axes;
axes('position', [0.15, 0.17, 0.8, 0.79]);
hold on;
set(gca, 'fontsize', 17); 
for m=1:scan_to_show
    tuning_curve = TC(:, m);
    if ~isempty(tuning_curve)
        plot(results.bins, tuning_curve, '-o', 'color', col_select(m,:));
    end
end
set(gca, 'xtick', 0:pi/2:2*pi);
set(gca, 'xticklabel', {'0', '\pi/2', '\pi', '3\pi/2', '2\pi'});
xlabel('Direction'); 
ylabel('Response'); 
xlim([0, 2*pi]);
box on;

if save_image
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_tuning_curve', em_id)));
    export_fig(gcf, fullfile(output_figs, sprintf('EM_%d_tuning_curve.pdf', em_id)));
end
