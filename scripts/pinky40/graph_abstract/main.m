%% load data 
ease.scan_id = 1; 
ease.block_id = 0;
ease.use_denoise = 0;
% load results 
neuron = ease.get_MF3D(false); 
% load CI data 
% Y = ease.load_Y(); 
% load EM 
% [Aem, segment_ids] = ease.load_Aem(); 

% prepare for exporting the results
save_image = true; 
output_figs = fullfile(ease.fig_folder, 'graph_abstract');
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end

%% example neurons 
cols = {'r', 'g'};  % colors for each neuron
em_ids = cell2mat(neuron.match_status.em_ids); 
em_id = 78665518; 
cell_idx = find(em_ids==em_id); 

%% show EM meshes  
ease.visualize_em_mesh(int64(em_id), true)
box on; 
set(gcf, 'color', 'w'); 
set(gca, 'boxstyle', 'full', 'xtick', [], 'ytick', [], 'ztick', [], ...
    'xlabel', [], 'ylabel', [], 'zlabel', []); 
axis equal tight; 

xvals = get(gca, 'xlim'); 
yvals = get(gca, 'ylim'); 
zvals = get(gca, 'zlim'); 
plot3([-25,-5]+xvals(2), [-10,-10]*0+yvals(1), [10, 10]+zvals(1), 'k', 'linewidth', 5);
if save_image
    export_fig(gcf, fullfile(output_figs, 'em_segment.fig'));
    export_fig(gcf, fullfile(output_figs, 'em_segment.pdf'));
    export_fig(gcf, fullfile(output_figs, 'em_segment.png'));
end

%% show tuning curves 
figure('papersize', [5, 3]); 
init_fig; 
tc = neuron.tuning_curve; 
plot(tc.x, tc.y(:, cell_idx), '-o', 'color', 'r', 'linewidth', 2, 'markersize', 10);
set(gca, 'position', [0.1, 0.1, 0.85, 0.85]); 
xlabel('Direction'); 
ylabel('Response'); 
axis off; 
if save_image
    export_fig(gcf, fullfile(output_figs, 'tuning_curve.fig'));
    export_fig(gcf, fullfile(output_figs, 'tuning_curve.pdf'));
    export_fig(gcf, fullfile(output_figs, 'tuning_curve.png'));
end

%% show temporal traces 
figure('papersize', [5, 1]); 
init_fig; 
T = size(neuron.C, 2); 
t = (1:T) / neuron.Fs;
plot(t, neuron.C_raw(cell_idx, :), 'color', 'b');
xlim([850, 900]);
ylim([-2, 70]); 
hold on; 
plot([895, 900], [15, 15], 'k', 'linewidth', 3); 
axis off; 
set(gca, 'position', [0.03, 0.1, 0.9, 0.9]); 
if save_image
    export_fig(gcf, fullfile(output_figs, 'temporal.fig'));
    export_fig(gcf, fullfile(output_figs, 'temporal.pdf'));
    export_fig(gcf, fullfile(output_figs, 'temporal.png'));
end
%% show spatial footrints
img1 = neuron.reshape(neuron.A(:, cell_idx), 3); 
[yy, xx] = determine_bounding_box(img1>max(img1(:))* 0.05, 0); 
figure('papersize', [diff(xx)+1, diff(yy)+1]/10); 
init_fig; 

imagesc((img1(yy(1):yy(2), xx(1):xx(2), 2)).^0.6); 
axis equal off; 
set(gca, 'position', [0, 0, 1, 1]); 
colormap jet; 
if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial.fig'));
    export_fig(gcf, fullfile(output_figs, 'spatial.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial.png'));
end

%% spatial mask EM
img2 = neuron.reshape(neuron.A_em(:, cell_idx), 3);
imwrite(img2(:, :, 2), fullfile(output_figs, 'em_mask.tif')); 

figure('papersize', [diff(xx)+1, diff(yy)+1]/10); 
init_fig; 

imagesc((img2(yy(1):yy(2), xx(1):xx(2), 2))); 
axis equal off; 
set(gca, 'position', [0, 0, 1, 1]); 
colormap jet; 
if save_image
    export_fig(gcf, fullfile(output_figs, 'spatial_em.fig'));
    export_fig(gcf, fullfile(output_figs, 'spatial_em.pdf'));
    export_fig(gcf, fullfile(output_figs, 'spatial_em.png'));
end
    
%% example frames 
ind_frames = 12113 + (1:3:30); 
dl = ease.video_loader{1, 2, 2}; 
x = dl.load_tzrc([3013, 3043]);
x = double(x(:, :, 1:3:end)); 
x = uint16(x / max(x(:)) * 65536*2); 
% writeTiff(uint16(x), fullfile(output_figs, 'example_frames.tif')); 
for m=1:size(x, 3)
    imwrite(x(:, :, m), fullfile(output_figs, sprintf('frame_%d.tif', m))); 
end 

%% example frames 
example_frames = Y(:, 12113+(1:3:30)); 
figure('papersize', [diff(xx)+1, diff(yy)+1]/10); 
init_fig;
for m=1:size(example_frames,2)
    img3 = neuron.reshape(example_frames(:, m), 3);
    %     imagesc(img3(yy(1):yy(2), xx(1):xx(2), 2), [0, 10000]);
    imagesc(img3);
    axis equal off;
    set(gca, 'position', [0, 0, 1, 1]);
    colormap gray;
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('frame_%d.fig', m)));
        export_fig(gcf, fullfile(output_figs, sprintf('frame_%d.pdf', m)));
        export_fig(gcf, fullfile(output_figs, sprintf('frame_%d.png', m)));
    end
end

%% create stimulus 
load example_stimuli.mat; 

[d1, d2, d3]  = size(img);

figure('papersize', [d2, d1]/100, 'position', [100, 100, d2, d1], 'unit', 'pixel'); 
imagesc(img); 
colormap gray; 
set(gca, 'position', [0, 0, 1, 1]); 
axis equal off tight; 

export_fig(gcf, fullfile(output_figs, 'stimulus.fig')); 
export_fig(gcf, fullfile(output_figs, 'stimulus.pdf'));
export_fig(gcf, fullfile(output_figs, 'stimulus.png')); 