%% order neurons based on spatial similarities 
neuron.orderROIs('temporal_cluster'); 
neuron.normalize('c_noise'); 

%%
K = size(neuron.A, 2); 
T_corr_spatial = cosine_similarity(neuron.A) - eye(K); 
T_corr_temporal = cosine_similarity(neuron.C') -eye(K); 

ind = true(size(T_corr_spatial)); 

% spatial correlation 
figure('papersize', [4.8, 4]); 
init_fig; 
axes('position',[0, 0.03, 0.99, 0.94]); 
imagesc(T_corr_spatial, [0, 1]); 
axis equal off tight; 
colorbar; 
colormap jet; 

if save_image
    export_fig(gcf, fullfile(output_figs, 'temporal_cluster_spatial.pdf'));
    export_fig(gcf, fullfile(output_figs, 'temporal_cluster_spatial.fig'));
else
    pause;
end

% temporal correlation 
figure('papersize', [4.8, 4]); 
init_fig; 
axes('position',[0, 0.03, 0.99, 0.94]); 
imagesc(T_corr_temporal, [0, 1]); 
axis equal off tight; 
colorbar; 
colormap jet; 

if save_image
    export_fig(gcf, fullfile(output_figs, 'temporal_cluster_temporal.pdf'));
    export_fig(gcf, fullfile(output_figs, 'temporal_cluster_temporal.fig'));
else
    pause;
end

%% known pairs 
% known_em_pairs = {[ 121738947, 96024719], ...
%     [92094831, 91299685], ...
%     [ 91693091, 53349510], ...
%     [ 75362044, 81004881], ...
%     [70018614,   55282126,   58281524], ...
%     [89057384, 77796334]}; 
known_em_pairs = [70018614,   55282126,   58281524];
npair = length(known_2p_pairs);
known_2p_pairs = 0*known_em_pairs;

em_ids = cell2mat(neuron.match_status.em_ids);
for n=1:length(known_em_pairs)
    tmp_idx = find(em_ids==known_em_pairs(n));
    known_2p_pairs(n) = tmp_idx;
end
C_ = neuron.C(known_2p_pairs, :);
C_raw = neuron.C_raw(known_2p_pairs, :);
c_merge = mean(C_, 1);

%% load the whole volume 
FOV = ease.FOV; 
ind_r = FOV(1):FOV(2); 
ind_c = FOV(3):FOV(4); 
dims_video = ease.dims_video; 
d1 = dims_video(1); 
d2 = dims_video(2); 
d3 = 3; 

if ~exist('corr_Yc.mat', 'file')
    Y_large = cell(d3, 1);
    neuron_large = MF3D('d1', d1, 'd2', d2, 'd3', d3);
    for mplane=1:3
        tmp_loader = ease.video_loader{ease.scan_id, mplane, ease.block_id};
        tmpY = tmp_loader.load_tzrc(neuron.frame_range);
        tmpY = reshape(double(tmpY), d1*d2, []);
        Y_large{mplane} = zscore(tmpY, 0, 2); %, d1*d2, []);
    end
    Y_large = cell2mat(Y_large);
        
    % compute correlation with the raw movie
    corr_Yc = Y_large*zscore(c_merge') / length(c_merge);
    save corr_Yc corr_Yc; 
else
    load corr_Yc.mat;
end

%% 
% extract neuron shapes and determine bounding box 
colors = [1, 0, 0; 0,1,1; 0, 1, 0];
A_ = imfilter(neuron.reshape(neuron.A(:, known_2p_pairs), 3), fspecial('gaussian', 3, 1)).^2; 
A_ = neuron.reshape(A_, 1); % enchance neural signals
for m=1:npair
   temp = A_(A_(:, m)>0, m); 
   A_(:, m) = A_(:,m) / quantile(temp(temp>0), 0.98) * 255; 
end
A_ = neuron.reshape(A_*colors(1:npair,:), 3); 
[rrange, crange] = determine_bounding_box(sum(A_, 4), 3); 
crop_r = round(mean(rrange)) + (-30:10) + FOV(1); 
crop_c = round(mean(crange)) + (-50:30) + FOV(3); 


% fuse them
img_bg = neuron_large.reshape(corr_Yc/max(corr_Yc(:)), 3); 
img = cell(3, 1);
for m=1:3 % slice by slice 
    % add neurons spatial footprints to the background
    temp = repmat(img_bg(:, :, m), [1, 1, 3]);
    temp = uint8(temp*255);
    
    % fuse
    for n=1:length(known_2p_pairs) % cell by cell 
        ai = squeeze(A_(:, :, m, :));        
        temp(ind_r, ind_c, :) = temp(ind_r, ind_c, :) + uint8(ai); 
    end 
    img{m} = 255-temp(crop_r, crop_c, :);
end

pixel_size = ease.range_2p(1) / ease.dims_video(1); 
neuron_large.showImage(img, 'horizontal', [], pixel_size, 'k'); 
tmp_axes = get(gcf, 'children'); 
em_ranges = ease.em_ranges(ease.video_zvals_updated(ease.scan_id,:)); 
for m=1:length(tmp_axes)
    axes(tmp_axes(m)); hold on;
    tmp_xlim = get(gca, 'xlim');
    tmp_ylim = get(gca, 'ylim');
    % draw contour
    xi = em_ranges{m}(:,2) / 2 ;
    yi = em_ranges{m}(:,1) /2 ;
    plot(xi+FOV(3)-crop_c(1), yi+FOV(1)-crop_r(1), '-.g', 'linewidth', 2);
    
    set(gca, 'xlim', tmp_xlim);
    set(gca, 'ylim', tmp_ylim);
end 
if save_image
     export_fig(gcf, fullfile(output_figs, 'example_merge.pdf'));
    export_fig(gcf, fullfile(output_figs, 'example_merge.fig'));
end 

%% save temporal traces 
n_merge = length(known_2p_pairs); 
y = neuron.C_raw(known_2p_pairs, :)/15;
y_denoised = neuron.C(known_2p_pairs, :)/15;

trange = 450+[0, 100];
tt = (1:size(y, 2)) / neuron.Fs;
ind = (tt>trange(1)) & (tt<=trange(2));

% temporal
figure('papersize', [6, 2]);
init_fig;
axes('position', [0.005, 0.005, 0.95, 0.99]); hold on;
set(gca, 'fontsize', 13);

for m=1:n_merge
    sn = std(y(m,:) - y_denoised(m,:), 0, 2); 
%     plot(tt(ind), y(m, ind)+ m, '-.', 'linewidth', 1, 'color', [1,1,1]*0.3 ); %color_list(m, :), 'linewidth', 1);
    plot(tt(ind), y(m, ind)+ m, 'color', 1-colors(m, :), 'linewidth', 2);
end
box on;
axis tight;
set(gca, 'xtick', []);
xlim(trange);
set(gca, 'fontsize', 18);
set(gca, 'xtick', []);
axis tight;
xlim(trange);
ylim([0.6, n_merge+1.3]);
set(gca, 'ytick', 1:n_merge); 
set(gca, 'yaxislocation', 'right'); 
box on;
% scale bar
plot(trange(end) + [-15, -5], [1, 1]*(0.73), 'k', 'linewidth', 4);

if save_image
    export_fig(gcf, fullfile(output_figs, 'example_merge_temporal.pdf'));
    export_fig(gcf, fullfile(output_figs, 'example_merge_temporal.fig'));
else
    pause;
end

%% show the complete neuron

% show mesh 
z_zoomin = 2.5;
% close all; 
% if ~exist('complete_mesh.mat', 'file')
%     [vertices, faces] = fetch1(ta3p100.Mesh & 'segmentation=2' & ...
%         'segment_id=648518346349483124', 'vertices', 'triangles');
%     load(fullfile(ease.dir_project, 'data/pinky100/coor_convert.mat'));
%     save complete_mesh.mat vertices faces A_convert offset;
% else
%     load complete_mesh.mat; 
% end
% 
% figure;
% vert = bsxfun(@plus, vertices * A_convert*0.001, offset);
% x = trisurf(faces+1, vert(:,1),...
%     ease.range_2p(2)-vert(:,2), ...
%     z_zoomin*(ease.range_2p(3)-vert(:,3)), ...
%     'edgecolor', 'none', 'facecolor', [0.25, 0.4, 0], 'facealpha', 0.3);
% hold on;

% add grid for specifying imaging planes 
zvals = ease.video_zvals_updated(ease.scan_id,:);
zvals = reshape(zvals', [], 1);
range_2p = ease.range_2p;
spatial_res = ease.range_2p(1:2)./ease.dims_video; 
spatial_range = neuron.reshape(neuron.spatial_range, 3);
ssub = 3; 
figure; hold on; 
for m = 1:3
    mask = double(spatial_range(:, :, m));
    mask(mask==0) =nan;
    [d1_crop, d2_crop] = size(mask);
    [xx, yy] = meshgrid(1:ssub:d2_crop, 1:ssub:d1_crop);
    yy = (yy+FOV(1)-1) * spatial_res(1);
    xx = (xx+FOV(3)-1) * spatial_res(2);
    surf(xx, yy, ones(size(xx))*zvals(m)*z_zoomin, mask(1:ssub:end, 1:ssub:end),...
        'edgecolor', [1, 1, 1]*0.8, 'edgealpha', 0.3, ...
        'facealpha', 0.3);
end
 
%show meshes of 3 divided segmented  
for m=1:length(known_em_pairs)
    ease.visualize_em_mesh(known_em_pairs(m), false, 1-colors(m,:), z_zoomin); 
    hold on; 
end
 
%
caxis([0,1]);
colormap([0, 1, 0]);
axis tight equal;
box on; 

xlim(crop_c([1,end])*spatial_res(2));
ylim(crop_r([1,end])*spatial_res(1));
zlim(z_zoomin*([min(zvals), max(zvals)]+[-12, 12]));
set(gca, 'xtick', 80:40:200);
set(gca, 'ytick', 160:40:300);
set(gca, 'ztick', (100:10:300)*z_zoomin, 'zticklabel', (100:10:300));
box on;
set(gca, 'boxstyle', 'full');
xlabel('X (um)', 'rotation', 7);
ylabel('Y (um)', 'rotation', -10);
zlabel('Z (um)');
set(gca, 'xdir', 'reverse'); 
set(gca, 'fontsize', 9); 
set(gcf, 'position', [300, 400, 500, 400]); 
set(gcf, 'papersize', [5, 4], 'color', 'w');
 set(gca, 'position', [0.13, 0.1, 0.85, 0.9]);
view(145, 8); 


if save_image
    print(gcf, fullfile(output_figs, 'example_merge_meshes.pdf'), '-dpdf', '-fillpage');
    export_fig(gcf, fullfile(output_figs, 'example_merge_meshes.png'));
    export_fig(gcf, fullfile(output_figs, 'example_merge_meshes.fig'));
else
    pause;
end
%% 
















