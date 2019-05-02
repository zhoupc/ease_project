%% order neurons based on spatial similarities
neuron.orderROIs('spatial_cluster');
neuron.normalize('c_noise');  %normalize data to make maximum of ci to be 1

example_em_ids = [
    128218233
    30242054
    63033382
    81004881];
K = size(neuron.A, 2); 
matched_em = cell2mat(neuron.match_status.em_ids);
example_idx = zeros(size(example_em_ids));
n_example = length(example_em_ids);
for m=1:length(example_em_ids)
    example_idx(m) = find(matched_em==example_em_ids(m));
end

%% if they are not together, reorder them 
if std(example_idx) ~= std(1:4)
    srt = 1:size(neuron.A, 2);
    new_idx = max(example_idx)+(1:4); % new id 
    srt(example_idx) = new_idx; 
    srt(new_idx) = example_idx;
    neuron.orderROIs(srt); 
    example_idx = new_idx; 
end
xmin = min(example_idx) - 0.5;
xmax = max(example_idx) + 0.5;
color_list = hsv(n_example);
n_example = length(example_idx);
d3 = neuron.options.d3;

%av%
K = size(neuron.A, 2);
S_corr_spatial = cosine_similarity(neuron.A) - eye(K);
S_corr_temporal = cosine_similarity(neuron.C') -eye(K);

ind = true(size(S_corr_spatial));

% spatial correlation
figure('papersize', [4.8, 4]);
init_fig;
axes('position',[0, 0.03, 0.94, 0.94]);
imagesc(S_corr_spatial, [0, 0.8]);
axis equal off tight;
colorbar;
colormap jet;

hold on;
plot([xmin, xmax, xmax, xmin, xmin], [xmin, xmin, xmax, xmax, xmin], 'r', 'linewidth', 2);
set(gca, 'fontsize', 15); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_spatial.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_spatial.fig'));
else
    pause;
end

% crop a small area
img_crop = S_corr_spatial(example_idx, example_idx);
figure('papersize', [2,2]);
init_fig;
axes('position', [0.18, 0.18, 0.78, 0.78]);
imagesc(img_crop, [0, 0.8]);
hold on;
plot([xmin, xmax, xmax, xmin, xmin]-xmin+0.5,...
    [xmin, xmin, xmax, xmax, xmin]-xmin+0.5, 'r', 'linewidth', 4);

colormap jet;
axis equal tight;
set(gca, 'xtick', [1, 2, 3, 4]);
set(gca, 'ytick', [1, 2, 3, 4]);
if save_image
    saveas(gcf, fullfile(output_figs, 'spatial_cluster_spatial_crop.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_spatial_crop.fig'));
else
    pause;
end

%% temporal correlation
figure('papersize', [4.8, 4]);
init_fig;
axes('position',[0, 0.03, 0.94, 0.94]);
imagesc(S_corr_temporal, [0, 0.8]);
axis equal off tight;
colorbar;
colormap jet;

hold on;
plot([xmin, xmax, xmax, xmin, xmin], [xmin, xmin, xmax, xmax, xmin], 'r', 'linewidth', 2);
set(gca, 'fontsize', 15); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_temporal.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_temporal.fig'));
else
    pause;
end

% crop a small area
img_crop = S_corr_temporal(example_idx, example_idx);
figure('papersize', [2,2]);
init_fig;
axes('position', [0.18, 0.18, 0.78, 0.78]);
imagesc(img_crop, [0, 0.8]);
hold on;
plot([xmin, xmax, xmax, xmin, xmin]-xmin+0.5,...
    [xmin, xmin, xmax, xmax, xmin]-xmin+0.5, 'r', 'linewidth', 4);

colormap jet;
axis equal tight;
set(gca, 'xtick', [1, 2, 3, 4]);
set(gca, 'ytick', [1, 2, 3, 4]);
if save_image
    saveas(gcf, fullfile(output_figs, 'spatial_cluster_temporal_crop.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_temporal_crop.fig'));
else
    pause;
end
%% determine spatial ranges
a = 5;
img = neuron.reshape(sum(neuron.A_em(:, example_idx), 2), 3);
[rrange, crange] = determine_bounding_box(sum(img, 3), 0);
temp = imfilter(sum(img,3), ones(a));
[r, c] = find(temp==max(temp(:)));
marker_x = c - crange(1);
marker_y = r - rrange(1);
w = (diff(crange) + 1) / 20;
h = diff(rrange) + 1 / 20;
theta = linspace(0, 2*pi, 80);
pixel_size = ease.range_2p(1) / ease.dims_video(1);
for m=1:n_example
    img_ai = neuron.reshape(neuron.A(:, example_idx(m)), 3);
    img_ai = img_ai(rrange(1):rrange(2), crange(1):crange(2), :);
    img_pi = neuron.reshape(neuron.A_em(:, example_idx(m)), 3);
    img_pi = img_pi(rrange(1):rrange(2), crange(1):crange(2), :);
    if m~=n_example
        neuron.showImage(img_ai, 'horizontal'); colormap jet;
    else
        neuron.showImage(img_ai, 'horizontal', [], pixel_size); colormap jet;
    end
    % plot markers
    temp = get(gcf, 'children');
    for n=1:length(temp)
        axes(temp(n));
        hold on;
        %         plot(marker_x, marker_y, '+w', 'markersize', 8);
        plot(marker_x+a*cos(theta), marker_y+a*sin(theta), '.w', 'linewidth', 2);
    end
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_2p_%d.pdf', m)));
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_2p_%d.fig', m)));
    else
        pause;
    end
    
    neuron.showImage(img_pi, 'horizontal', []); colormap jet;
    % plot markers
    temp = get(gcf, 'children');
    for n=1:length(temp)
        axes(temp(n));
        hold on;
        %         plot(marker_x, marker_y, '+w', 'markersize', 8);
        plot(marker_x+a*cos(theta), marker_y+a*sin(theta), '.w', 'linewidth', 2);
    end
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_em_%d.pdf', m)));
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_em_%d.fig', m)));
    else
        pause;
    end
end

%%
y = neuron.C_raw(example_idx, :)/15;
y_denoised = neuron.C(example_idx, :)/15;
trange = [50, 150];
tt = (1:size(y, 2)) / neuron.Fs;
ind = (tt>trange(1)) & (tt<=trange(2));

% temporal
% temporal
figure('papersize', [6.5, 2]);
init_fig;
axes('position', [0.005, 0.005, 0.97, 0.99]); hold on;
set(gca, 'fontsize', 13);

for m=1:n_example
    sn = std(y(m,:) - y_denoised(m,:), 0, 2); 
    plot(tt(ind), y(m, ind)+ m, '-.', 'linewidth', 1, 'color', [1,1,1]*0.3 ); %color_list(m, :), 'linewidth', 1);
    plot(tt(ind), y_denoised(m, ind)+ m, 'color', color_list(m, :), 'linewidth', 2);
end
box on;
axis tight;
set(gca, 'xtick', []);
xlim(trange);
set(gca, 'fontsize', 14);
set(gca, 'xtick', []);
axis tight;
xlim(trange);
ylim([0.6, 5]);
set(gca, 'ytick', 1:4); 
set(gca, 'yaxislocation', 'right'); 
box on;
% scale bar
plot(trange(end) + [-15, -5], [1, 1]*(0.73), 'k', 'linewidth', 4);


%
% axes('position', [0.1, 0, 0.88, 0.1]);
% plot([175, 195], [0.7, 0.7], 'k', 'linewidth', 8);
% xlim([0, 200]);
% axis  off;
% ylim([0, 1]);
if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_example_temporal.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_cluster_example_temporal.fig'));
else
    pause;
end

%% create a video 
Y_res = neuron.compute_residual(Y); 
sn_pixels = std(Y_res, 0, 2); 
Y_res = bsxfun(@times, Y_res, 1./sn_pixels); 
Y_res(~neuron.spatial_range, :) = nan; 

A_show = bsxfun(@times, full(neuron.A(:, example_idx)), 1./sn_pixels); 
C_show = neuron.C(example_idx, :); 
S_show = neuron.S(example_idx, :); 
A_em = full(neuron.A_em(:, example_idx)); 

Y_show = Y_res + A_show*C_show; 

%%
% rotate the FOV 
[rot_angle, rot_xlim, rot_ylim] = find_rotation(neuron.reshape(sum(neuron.A_em, 2),3)); 
A_em = neuron.reshape(A_em, 3); 
[d1, d2, d3, n_example] = size(A_em); 
A_em_rot = imrotate(reshape(A_em, d1, d2, []), rot_angle); 
[rrange, crange] = determine_bounding_box(A_em_rot, 0); 
rrange = rrange + [2, -2]; 
crange = crange + [20, -3]; 
d1_rot = diff(rrange) + 1; 
d2_rot = diff(crange) + 1; 
A_em_rot = A_em_rot(rrange(1):rrange(2), crange(1):crange(2), :); 
A_em_rot = reshape(A_em_rot, d1_rot, d2_rot, d3, []); 

A_show_rot = imrotate(reshape(A_show, d1, d2, []), rot_angle); 
A_show_rot = A_show_rot(rrange(1):rrange(2), crange(1):crange(2), :); 
A_show_rot = reshape(A_show_rot, d1_rot, d2_rot, d3, []); 

Y_show_rot = imrotate(reshape(Y_show, d1, d2, []), rot_angle); 
Y_show_rot = Y_show_rot(rrange(1):rrange(2), crange(1):crange(2), :); 
Y_show_rot = reshape(Y_show_rot, d1_rot, d2_rot, d3, []); 

Y_res_rot = imrotate(reshape(Y_res, d1, d2, []), rot_angle); 
Y_res_rot = Y_res_rot(rrange(1):rrange(2), crange(1):crange(2), :); 
Y_res_rot = reshape(Y_res_rot, d1_rot, d2_rot, d3, []); 

%% compute contours for all 4 neurons 

cell_contours = cell(n_example, 1); 
for m=1:n_example 
    cell_contours{m} = get_contour(A_show_rot(:, :, :, m), 0.75); 
end 

height = 727; 
width = 554; 
figure('position', [512, 412, width, height]);
ax = tight_subplot(8, 3); 
ax = reshape(ax, 3, 8)'; 
h_img = cell(8, 3); 
% merge the last 6 panels into one 
temp = ax(8, 1).Position; 
ax_pos = temp; 
temp = ax(8, 3).Position; 
ax_pos(3) = temp(1)+temp(3) - ax_pos(1); 
temp = ax(7, 1).Position; 
ax_pos(4) = temp(2) + temp(4) - ax_pos(2);
ax_pos(1) = ax_pos(1) *1.5; 
ax_temporal = axes('position', ax_pos); 
hold(ax_temporal, 'on'); 
ylim_temporal = [0.77, 5]; 
set(ax_temporal, 'ylim', ylim_temporal); 

% draw contours
temp = sum(sum(A_em_rot, 3),4);
[r, c] = find(temp==max(temp(:)));
marker_x = c;
marker_y = r;
theta = linspace(0, 2*pi, 80);

for m=1:8
    for n=1:3
        
        if m>6
            delete(ax(m, n));
            continue;
        end
        h_img{m, n} = imagesc(temp, 'parent', ax(m, n));
        hold(ax(m,n), 'on');
        plot(marker_x+a*cos(theta), marker_y+a*sin(theta), ...
            'y', 'linewidth', 0.5, 'parent', ax(m, n));
        axis(ax(m,n), 'off', 'tight', 'equal');
    end
end
for m=1:n_example
    for n=1:3
        cont = cell_contours{m}{n};
        for k=1:length(cont)
            temp = cont{k}; 
           plot(temp(1,:), temp(2,:), 'color', color_list(m,:), ...
                'linewidth', 1,'parent', ax(m+2, n)); 
        end
    end  
end

color_text = 'm'; 
text(1, -2, 'subtracted raw signal', 'parent', ax(1,1), ...
    'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
text(1, -2, 'residual', 'parent', ax(2,1), 'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
text(1, -2, 'component 1', 'parent', ax(3,1), 'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
text(1, -2, 'component 2', 'parent', ax(4,1), 'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
text(1, -2, 'component 3', 'parent', ax(5,1), 'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
text(1, -2, 'component 4', 'parent', ax(6,1), 'color', color_text, 'fontweight', 'bold', 'fontsize', 10); 
colormap bone; 
t_now = 0; 
tmp_h = text(1, d1_rot-2, sprintf('Time = %.2f sec', t_now),...
    'fontsize', 10, 'color', 'y', 'parent', ax(6,1));
for m=1:n_example
    sn = std(y(m,:) - y_denoised(m,:), 0, 2); 
    plot(tt, y(m,:)+ m, '-.', 'linewidth', 1, 'color', [1,1,1]*0.3, ...
        'parent', ax_temporal); 
    plot(tt, y_denoised(m, :)+ m, 'color', color_list(m, :),...
        'linewidth', 2, 'parent', ax_temporal);
end
time_line = plot([1,1]*10, ylim_temporal, '-.m',...
    'linewidth', 1,'parent', ax_temporal); 
    
box(ax_temporal, 'on'); 
xlabel(ax_temporal, 'Time (Sec.)'); 
ylabel(ax_temporal, 'Component #'); 
set(ax_temporal, 'xlim', [0, 50]); 

%% play video 
kt = 3;
avi_nm = fullfile(ease.video_folder, 'demixing_example.avi'); 
avi_file = VideoWriter(avi_nm);
if ~isnan(neuron.Fs)
    avi_file.FrameRate = neuron.Fs;
end
avi_file.Quality = 100;
avi_file.open();
avi_flag = true;

ind_active = (sum(C_show, 1)>0); 
for frame_idx = 1:kt:8900 % 1378, %7555 %1377
    img_show = Y_show_rot(:, :, :, frame_idx);
    img_res = Y_res_rot(:, :, :, frame_idx);
    for m=1:3
        % raw data
        delete(h_img{1,m});
        h_img{1, m} = imagesc(img_show(:, :, m),'parent', ax(1, m), [-4, 4]);
        temp = (get(ax(1,m), 'children'));
        set(ax(1,m), 'children', temp([2:end, 1]));
        
        % residual
        delete(h_img{2,m});
        h_img{2, m} = imagesc(img_res(:, :, m), 'parent', ax(2, m), [-4,4]);
        temp = (get(ax(2,m), 'children'));
        set(ax(2,m), 'children', temp([2:end,1]));
        
        % neurons
        for k=1:n_example
            delete(h_img{k+2,m});
            h_img{k+2, m} = imagesc(A_show_rot(:, :, m, k)*...
                C_show(k, frame_idx),'parent', ax(k+2, m), [-4,4]);
            temp = (get(ax(k+2,m), 'children'));
            set(ax(k+2,m), 'children', temp([2:end,1]));
        end
    end
    
    t_now = frame_idx/neuron.Fs; 
    delete(tmp_h); 
    tmp_h = text(1, d1_rot-2, sprintf('Time = %.2f sec', t_now),...
        'fontsize', 10, 'color', 'y', 'parent', ax(6,1));
    
    set(time_line, 'xdata', t_now*[1,1]);
    set(ax_temporal, 'xlim', max(0, t_now-25)+[0, 50]); 
    drawnow(); 
    
    if avi_flag
        temp = getframe(gcf);
        temp.cdata = imresize(temp.cdata, [height, width]);
        avi_file.writeVideo(temp);
    end
    delete(tmp_h); 
end
avi_file.close(); 










