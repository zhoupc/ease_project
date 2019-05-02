%% order neurons based on PNR
neuron.orderROIs('snr');
K = size(neuron.A, 2);
n_example = 20;
n_per_row = 4; 
ind_example = (1:n_example)*5;

%% extract neuron shapes
A_ = neuron.A(:, ind_example);
A_ = bsxfun(@times, A_, 1./max(A_, [],1));
A_em = neuron.A_em(:, ind_example);
A_em = bsxfun(@times, A_em, 1./max(A_em, [], 1));

% rotate & zoomin 
[~, rot_info] = neuron.rotate_zoomin(neuron.spatial_range); 
A_ = neuron.rotate_zoomin(A_, rot_info);
A_em = neuron.rotate_zoomin(A_em, rot_info);

%%
img_2p = squeeze(max(A_, [], 3));
img_em = squeeze(max(A_em, [], 3));
% img_2p = squeeze(sum(A_, 3));
% img_em = squeeze(sum(A_em, 3));

temp = size(img_2p);
x0 = temp(2)-18;
y0 = 6;
for m=1:n_per_row:(n_example-n_per_row+1)
    if m==n_example - n_per_row+1
        pixel_size = ease.range_2p(1) / ease.dims_video(1); 
    else
        pixel_size = []; 
    end
        
    %% show 2p neuron
    neuron.showImage(img_2p(:, :, m:(m+3)), 'horizontal', [], pixel_size);
    colormap jet;
    
%     tmp_ax = get(gcf, 'children');
%     for n=1:5
%         text(x0, y0, num2str(ind_example(m+n-1)), 'color', 'w', 'parent', tmp_ax(n),...
%             'fontweight', 'bold', 'fontsize', 12);
%     end
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('2p_example_%d_%d.pdf',m, m+n_per_row-1)));
        export_fig(gcf, fullfile(output_figs, sprintf('2p_example_%d_%d.fig',m, m+n_per_row-1)));
    else
        pause;
    end
    close; 
    
    %% show em neurons
    neuron.showImage(img_em(:, :, m:(m+3)),'horizontal');
    colormap jet;
    
    tmp_ax = get(gcf, 'children');
    for n=1:n_per_row
        text(x0, y0, num2str(ind_example(m+n-1)), 'color', 'w', 'parent', tmp_ax(n),...
            'fontweight', 'bold', 'fontsize', 11);
    end
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('em_example_%d_%d.pdf',m, m+n_per_row-1)));
        export_fig(gcf, fullfile(output_figs, sprintf('em_example_%d_%d.fig',m, m+n_per_row-1)));
    else
        pause;
    end
    close; 
end

%%
C_raw = full(neuron.C_raw(ind_example, :));
C_ = full(neuron.C(ind_example, :));
C_raw_scale= bsxfun(@times, C_raw, 0.6./std(C_raw-C_, 0, 2));
C_scale= bsxfun(@times, C_, 0.6./std(C_raw-C_, 0, 2));
color_list = cool(n_example);
% color_list = flipud(color_list);

figure('papersize', [8, 2.9]);
init_fig;
axes('position', [0.035, 0.01, 0.959, 0.985]);
hold on; 
set(gca, 'yaxislocation', 'left'); 
T = size(C_raw_scale, 2);
tt = (1:T) / neuron.Fs;
xmin = 50; 
xmax = 150; 
ymin = 0; 
ymax = ind_example(end)+5; 
ind_t = (tt>=xmin) & (tt<=xmax); 
for m=n_example:-1:1
%       plot(tt, C_raw_scale(m,:)+ ind_example(m), '-.', 'linewidth', 1, 'color', [1,1,1]*0.3 ); %color_list(m, :), 'linewidth', 1);
%     plot(tt, C_scale(m, :)+ ind_example(m), 'color', color_list(m, :), 'linewidth', 2);
% 
    plot(tt(ind_t), ind_example(m)+C_raw_scale(m,ind_t), '-.', 'color', [1,1,1]*0.5, 'linewidth', 2); %'color', color_list(m,:), 'linewidth', 0.5);
    plot(tt(ind_t), ind_example(m)+C_scale(m,ind_t), 'color', color_list(m,:), 'linewidth', 1);
ymax = max(ymax, ind_example(m)+max(C_raw_scale(m,ind_t))); 

end

set(gca, 'xtick', []);
set(gca, 'ytick', ind_example(n_per_row:n_per_row:n_example));
box off;
xlim([xmin, xmax]); 
ylim([ymin, ymax]);
box on; 
% plot([xmax, xmax, xmin, xmin, xmax], [ymin, ymax, ymax, ymin, ymin], 'k'); 
set(gca, 'fontsize', 10); 

% axes('position', [0.05, 0., 0.94, 0.05]);
plot([xmax-15, xmax-5], [1,1]*1.8, 'color', 'k', 'linewidth', 4);

if save_image
    export_fig(gcf, fullfile(output_figs, 'temporal_examples.pdf'));
    export_fig(gcf, fullfile(output_figs, 'temporal_examples.pdf'));
else
    pause;
end



