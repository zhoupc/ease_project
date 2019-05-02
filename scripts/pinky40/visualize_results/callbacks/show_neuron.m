cell_id = max(1,cell_id); 
cell_id = min(K, cell_id); 
set(input_cellid, 'string', num2str(cell_id)); 
%% show spatial footprints
id_2p = idx_show(cell_id);
ai = reshape(full(results.A(:, id_2p)), results.dims);
ai_corr = reshape(full(results.A_corr(:, id_2p)), results.dims);
ai_em = reshape(full(results.A_em(:, id_2p)), results.dims);
confidence = results.confidences(:, id_2p); 
labels = results.labels(:, id_2p); 
r2 = 1-results.rss(:, id_2p); 

for m=1:nscan
    
    for n=1:d3
        % show ai_corr
%         axes(ax_corr{n, m});
        imagesc(squeeze(ai_corr(:, :, n, m)), 'parent', ax_corr{n, m});
        axis(ax_corr{n, m}, 'off');
        
        % show ai
        imagesc(squeeze(ai(:, :, n, m)), 'parent', ax_ai{n, m});
        axis(ax_ai{n, m}, 'off');
        
        % show pi
        imagesc(squeeze(ai_em(:, :, n, m)), 'parent', ax_pi{n, m});
                axis(ax_pi{n, m}, 'off');

    end
    
    set(text_confidence{m}, 'string', sprintf('%d: %.3f', scan_list(m), confidence(m))); 
    if labels(m)<0
        set(text_confidence{m}, 'foregroundcolor', 'w', 'fontweight', 'bold', 'fontsize', 18);
    else
        set(text_confidence{m}, 'foregroundcolor', 'k', 'fontweight', 'normal', 'fontsize', 16);
    end
end

%% show temporal traces 
C_raw = results.Craw(:, :, id_2p); 
C = results.C(:, :, id_2p); 
cla(ax_temporal); 
hold(ax_temporal, 'on'); 
dy = 20; 
for m=1:nscan
    plot(C_raw(:,m)+dy*m, 'color', [1,1,1]*0.9, 'linewidth', 0.5, 'parent', ax_temporal); 
    hold on; 
end
for m=1:nscan
    plot(C(:,m)+dy*m, 'color', colors(m,:), 'linewidth', 2, 'parent', ax_temporal); 
end
set(gca, 'ytick', dy:dy:(nscan*dy)); 
set(gca, 'yticklabel', scan_list); 
axis(ax_temporal, 'tight'); 
set(ax_temporal, 'xlim', [0, 5000]); 

%% show tuning curves 
tc_y = results.y(:, :, id_2p); 
tc_yfit = results.yfit(:, :, id_2p); 
cla(ax_tc); 
hold(ax_tc, 'on'); 
for m=1:nscan 
    plot(results.bins, tc_y(:, m), '-o', 'color', colors(m,:), 'linewidth', 2); 
    set(text_r2{m}, 'string', sprintf('%.3f', r2(m))); 
    if labels(m)<0
        set(text_r2{m}, 'foregroundcolor', 'w', 'fontweight', 'bold', 'fontsize', 18);
    else
        set(text_r2{m}, 'foregroundcolor', 'k', 'fontweight', 'normal', 'fontsize', 16);
    end
end 
