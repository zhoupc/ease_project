%% prepare for running EASE 
if ~exist('ease', 'var')
    fi.usepkg('ease'); 
    run_ease; 
end

%% prepare for exporting the results 
save_image = true;      % export figures or not 

output_figs = fullfile(ease.fig_folder, 'EM_features');% folder for exporting figures 
if ~exist(output_figs, 'dir')
    mkdir(output_figs);
end

%% load EM footprints 
ease.scan_id = 1; 
neuron = ease.get_MF3D(true); 
[Aem, segment_ids] = ease.load_Aem(); 

%% example cells to show  
% example_ids =  uint64([90292018, ...
%     110587899, ...
%     80994419, ...
%     44523552, ...
%     110497688]); 
example_ids =  uint64([90292018, ...
    77796334, ...
    33376321, ...
    84930953, ...
    98023685]);

%% show histogram of Aem
v = sum(Aem, 1);
[v, ind_sort] = sort(v, 'descend');
v_max = max(v); 
K = length(v);

% find the rank of those example neurons 
example_rank = zeros(size(example_ids));
tmp_ids = segment_ids(ind_sort); 
for m=1:length(example_ids)
    example_rank(m) = find(tmp_ids==example_ids(m), 1, 'first'); 
end 

% figure & axes 
figure('papersize', [8, 3.5]);
init_fig;
axes('position', [0.08, 0.22, 0.91, 0.73]);

% plot the normalized values 
plot(1:K, v/v_max, 'k');
set(gca, 'yscale', 'log');
set(gca, 'xscale', 'log');
xlim([0, 45000]);
ylim([floor(v(45000)/10)*10, v(1)]/v_max);
hold on;
box on;

% add examples 
x_text = [0.8, 0.7, 0.7, 0.7, 0.22];  % adjust the position of labels 
y_text = [0.5, 2.2, 0.38, 2.6, 1];
markers = {'s', 'o', 'v', 'd', '*'};
col = {'r', 'g', 'c', 'm', 'b'};
for m=1:length(example_rank)
    k = example_rank(m);
    plot(k, v(k)/v_max, markers{m}, 'color', col{m}, 'markersize', 15);
    text(k*x_text(m), v(k) *y_text(m)/v_max, sprintf('%d', k), 'fontsize',...
        20, 'color', col{m});
end
xlabel('rank');
set(gca, 'xtick', [1, 10, 100, 1000, 10000]);
set(gca, 'xticklabel', [1, 10, 100, 1000, 10000]);
set(gca, 'ytick', [1000, 10000, 100000,1000000]/10^6);
if save_image
    export_fig(gcf, fullfile(output_figs, 'nvoxels.pdf'));
    export_fig(gcf, fullfile(output_figs, 'nvoxels.fig'));
end

% save spatial footprints 
for m=1:length(example_rank)
    k = example_rank(m);
    ai = neuron.reshape(Aem(:, ind_sort(k)), 3);
    
    if m==length(example_rank)
        pixel_size = ease.range_2p(1) / ease.dims_video(1); 
    else
        pixel_size = []; 
    end
    neuron.showImage(ai, 'horizontal', [], pixel_size);
     
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('example_%d.pdf', m)));
        export_fig(gcf, fullfile(output_figs, sprintf('example_%d.fig', m)));
    end
end







