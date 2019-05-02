white_bg = false; 
use_EM = false; 
thr_min = 0.1; 
spatial_range = neuron.reshape(neuron.spatial_range, 3);
[rot_angle, rot_xlim, rot_ylim] = find_rotation(spatial_range);
orientation = 'horizontal'; 

%% all somas
neuron.orderROIs('confidence'); 
K = size(neuron.A, 2); 
n_batch = 5; 
ind = [0, 10, 30, 60, 100, K];
k0 = ind(1)+1; 
for mbatch=1:n_batch
    k1 = ind(mbatch+1); 
    img = neuron.overlapA(k0:k1, thr_min, white_bg, use_EM);
    for m=1:3
        temp = imrotate(img{m}, rot_angle);
        img{m} = temp(rot_ylim(1):rot_ylim(2), rot_xlim(1):rot_xlim(2), :);
    end
    if mbatch==n_batch 
        use_scale_bar = pixel_size; 
    else
        use_scale_bar = []; 
    end
    neuron.showImage(img, orientation, [], use_scale_bar);
    if mbatch==n_batch 
        k1 = inf; 
    end
    if save_image
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_overlap_%d_%d_%d.pdf', k0, k1, white_bg)));
        export_fig(gcf, fullfile(output_figs, sprintf('spatial_overlap_%d_%d_%d.fig', k0, k1, white_bg)));
    end
    k0 = k1 + 1; 
end

disp(neuron.match_status.confidence(ind(2:end))); 
% 2.9061    2.1531    1.5719    1.2298    0.4028
