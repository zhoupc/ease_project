%% load results
addpath('callbacks');
% load results.mat ;
dims = results.dims;
scan_list = results.scan_list;
[nscan, K] = size(results.confidences);
max_confidence = max(results.confidences, [], 1);
[~, idx_show] = sort(max_confidence, 'descend');
cell_id = 1;
colors = cool(4); 

%%
close all;
w1 = dims(2)*0.8;
h1 = dims(1)*0.8;
dw1 = 10;
dw2 = 2;
w = 2 * dims(3)*(w1+dw1+dw2*(dims(3)-1));
h = 940;

w0 = 30;
h0 = h-h1-20;

% create a main figure
fig_main = figure('position', [10, 10, w, h], 'KeyPressFcn', @keyPress);

% select box
w2 = 40;
pos_cell_id_minus = [w0, h-20, w2, 20];
pos_cell_id = pos_cell_id_minus + [w2, 0, 0, 0];
pos_cell_id_plus = pos_cell_id + [w2, 0, 0, 0];
input_cellid = uicontrol('parent', fig_main, 'style', 'edit', ...
    'unit', 'pixel', 'position', pos_cell_id, 'fontsize', 16,...
    'string', num2str(cell_id), 'callback', ...
    'cell_id=str2double(get(gco, ''string'')); show_neuron');
cell_id_minus = uicontrol('parent', fig_main, 'style', 'pushbutton', ...
    'unit', 'pixel', 'position', pos_cell_id_minus, 'fontsize', 16,...
    'string', '<<', 'callback', ...
    'cell_id=max(1, cell_id-1); show_neuron');
cell_id_plus = uicontrol('parent', fig_main, 'style', 'pushbutton', ...
    'unit', 'pixel', 'position', pos_cell_id_plus, 'fontsize', 16,...
    'string', '>>', 'callback', ...
    'cell_id=min(K, cell_id+1); show_neuron');

% axes for showing spatial components
d3 = dims(3);
ax_corr = cell(d3, nscan);
ax_ai = cell(d3, nscan);
ax_pi = cell(d3, nscan);
for m=1:4
    tmp_w0 = w0+(1-mod(m, 2))*d3*(w1+dw1);
    tmp_h0 = h0 - (floor((m+1)/2)-1) * 3* (h1+10)-10;
    for n=1:d3
        tmp_pos = [tmp_w0+w1*(n-1)+dw2, tmp_h0, w1, h1];
        ax_corr{n, m}= axes('parent', fig_main, ...
            'unit', 'pixel', 'position', tmp_pos);
        
        tmp_pos(2) = tmp_pos(2) - h1-10;
        ax_ai{n, m}= axes('parent', fig_main, ...
            'unit', 'pixel', 'position', tmp_pos);
        
        tmp_pos(2) = tmp_pos(2) - h1-10;
        ax_pi{n, m}= axes('parent', fig_main, ...
            'unit', 'pixel', 'position', tmp_pos);
    end
end

% show temporal
temp = get(ax_pi{1, 3}, 'position');
pos_temporal = [temp(1), temp(2)-320, w-w0-10, 300];
ax_temporal = axes('parent', fig_main,...
    'unit', 'pixel', 'position', pos_temporal, 'fontsize', 16);

%% axes for visualizing tuning curves
pos_tc = [pos_temporal(1), 10, 410, pos_temporal(2)-40];
ax_tc = axes('parent', fig_main,...
    'unit', 'pixel', 'position', pos_tc, 'fontsize', 16);

%%
pos_confidence = [sum(pos_tc([1,3]))+10, sum(pos_tc([2, 4]))-20, 100, 20];
btn_confidence = uicontrol('parent', fig_main, 'style', 'pushbutton', ...
    'unit', 'pixel', 'position', pos_confidence, 'string', 'confidence', ...
    'fontsize', 16);

text_confidence = cell(nscan, 1);
for m=1:nscan
    text_confidence{m} = uicontrol('parent', fig_main, 'style', 'text', ...
        'unit', 'pixel', 'position', pos_confidence-[0, 30*m, 0, 0], ...
        'backgroundcolor', colors(m,:), 'fontsize', 16);
end
%

pos_r2 = pos_confidence + [pos_confidence(3)+10, 0, 0, 0];
btn_r2 = uicontrol('parent', fig_main, 'style', 'pushbutton', ...
    'unit', 'pixel', 'position', pos_r2, 'string', 'R-squared', ...
    'fontsize', 16);

text_r2 = cell(nscan, 1);
for m=1:nscan
    text_r2{m} = uicontrol('parent', fig_main, 'style', 'text', ...
        'unit', 'pixel', 'position', pos_r2-[0, 30*m, 0, 0], ...
        'backgroundcolor', colors(m,:), 'fontsize', 16);
end
% pos_confidence = [sum(pos_tc([1,3]))+10, sum(pos_tc([2, 4]))-20, 100, 20];
% btn_confidence = uicontrol('parent', fig_main, 'style', 'pushbutton', ...
%     'unit', 'pixel', 'position', pos_confidence, 'string', 'confidence', ...
%     'fontsize', 16);
% text_1 = uicontrol('parent', fig_main, 'style', 'text', ...
%     'unit', 'pixel', 'position', pos_confidence-[0, 30, 0, 0], ...
%     'backgroundcolor', colors(1,:), 'fontsize', 16, 'string', '1: 0');
% text_2 = uicontrol('parent', fig_main, 'style', 'text', ...
%     'unit', 'pixel', 'position', pos_confidence-[0, 60, 0, 0], ...
%     'backgroundcolor', colors(2,:), 'fontsize', 16, 'string', '2: 0');
% text_3 = uicontrol('parent', fig_main, 'style', 'text', ...
%     'unit', 'pixel', 'position', pos_confidence-[0, 90, 0, 0], ...
%     'backgroundcolor', colors(3,:), 'fontsize', 16, 'string', '3: 0');
% text_4 = uicontrol('parent', fig_main, 'style', 'text', ...
%     'unit', 'pixel', 'position', pos_confidence-[0, 120, 0, 0], ...
%     'backgroundcolor', colors(4,:), 'fontsize', 16, 'string', '4: 0');

%%
show_neuron;

