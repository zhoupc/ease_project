img = neuron.reshape(img, 3);
orientation = 'horizontal';
d1 = diff(yrange)+1;
d2 = diff(xrange)+1;

if strcmpi(orientation, 'horizontal')
    h = d1+2;
    w = d3*(d2+1);
    pos_ax = [1-(d2+1)/w, 1/h, d2/w, 1-2/h];
    dpos = [-(d2+1)/w, 0, 0, 0];
else
    h = d3*(d1+1);
    w = d2+1;
    pos_ax = [1/w, 1/h, 1-2/w, d1/h];
    dpos = [0, (d1+1)/h, 0, 0];
end
figure('papersize', [w, h] / max(w, h)*10);
init_fig;
mscan = ease.scan_id;
ssub = 2;
for m=1:3
    axes('position', pos_ax);
    imagesc(img(:, :, d3-m+1), vlim);
    axis off; %equal off tight;
    pos_ax = pos_ax + dpos;
    hold on;
    
    % draw contour
    if ~isempty(em_ranges{m})
        xi = em_ranges{m}(:,1);
        yi = em_ranges{m}(:,2);
        plot(xi, yi, '-.k', 'linewidth', 3);
    end
    xlim(xrange);
    ylim(yrange);
    
    % add a scale bar
    if exist('add_scale_bar', 'var') && add_scale_bar && (m==1)
        plot(xrange(1)+d2-3- [0, 20]/pixel_size, ...
            yrange(1)+d1-[1,1]*3, 'w', 'linewidth', 5);
    end
end
colormap jet;

