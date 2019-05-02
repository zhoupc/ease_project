disp('you should run this script with GUI.'); 

%% fix an ambguous match. in practice, we don't do this step, but here we need to reproduce the paper results 
em_ids = cell2mat(neuron.match_status.em_ids); 
ind = find(em_ids== 124640933);
if ~isemtpy(ind)
    neuron.match_status.em_ids{ind} = 128218233; 
end 
%% delete neurons that stealed signals from neighboring neurons
neuron.orderROIs('spatial_cluster'); 
correlated_pairs = neuron.find_correlated_pairs([0.1, 0.3]);
ease.startGUI(); 

for m=1:size(correlated_pairs, 1)
    tmp_pair = correlated_pairs(m,:); 
    neuron.compare_pairs(tmp_pair); 
    disp(tmp_pair);
    fprintf('which one to delete? (green: %d, magenta: %d)\n', tmp_pair(1), tmp_pair(2)); 
    temp = input('choose (g: green; m: magenta; n: none, b: both): ', 's');
    if strcmpi(temp, 'g')
        neuron.match_status.confidence(tmp_pair(1)) = 0;
    elseif strcmpi(temp, 'm')
        neuron.match_status.confidence(tmp_pair(2)) = 0;
    elseif strcmpi(temp, 'b')
        neuron.match_status.confidence(tmp_pair) = 0;
    end
    ease.closeall(); 
end
neuron.delete(neuron.match_status.confidence==0, true)
neuron.update_spatial(Y); 
% neuron.update_temporal(Y); 

%% delete bad quality neurons 
neuron.evaluate_matching_confidence(Y, Aem, segment_ids); 
neuron.orderROIs('confidence'); 
ease.startGUI(); 
pause; 
neuron.merge_repeats();
neuron.delete(neuron.match_status.confidence==0, true); 
neuron.hals(Y); 

%% to reproduce the results in the paper, we create a white list of neurons. usually these neurons were initialized in the earlier steps 
white_list = [
    70018614  % example for showing temporally correlated segments 
    55282126
    58281524     
    128218233
    30242054
    63033382
    81004881
    90292018  % somas 
    88764136
    23410728
    79594839
    83286327
    82868719
    94050687
    82251040
    78665518
    50713226
    50725796
    71868389
    73839508
    59520759
   101666914
    81578447
    36166162
    93217930
    73430227
    92481124
    52089750
    84832010
    88240816
    54548179
    77796334
    ]; 
options_init = ease.options_init; 
options_init.clear_results = false; 
options_init.show_fig = false; 
options_init.K_candidate = 40000; 
neuron.ease_initialization(Y, Aem, segment_ids, options_init, [], white_list);
neuron.hals(Y); 
neuron.evaluate_matching_confidence(Y, Aem, segment_ids); 

%% label neurons 
K = size(neuron.A, 2); 
labeled_cell = [90292018, -1
79594839, -1
10914339, -2
20410486, -2
23410728, -1
88764136, -1
76859794, -2
78665518, -1
83286327, -1
90703146, -1
81004881, -2
82868719, -1
18572815, -2
121022936, -2
101666914, -1
36166162, -1
17660351, -2
26113030, -2
77507789, -2
52089750, -1
81578447, -1
71868389, -1
19315505, -2
82251040, -1
59520759, -1
114090833, -2
67981253, -2
121738947, -2
50713226, -2
84832010, -1
94050687, -1
55282126, -2
5386310, -2
93217930, -1
73430227, -1
73839508, -1
30242054, -2
54548179, -1
115217795, -2
82848773, -2
107860456, -2
91693091, -2
91299685, -2
77796334, -1
51320039, -2
88240816, -1
92481124, -1
63033382, -2
50725796, -1
99046706, -2
28637836, -2
694582, 2
5280949, -2
92948535, -2
107731759, -2
48234238, -2
58605699, -2
58281524, -2
60064737, -2
66936105, -2
70018614, -2]; 

em_ids = cell2mat(neuron.match_status.em_ids);
for m=1:size(labeled_cell,1)
    k = find(em_ids==labeled_cell(m, 1));
    neuron.labels(k)=labeled_cell(m,2);
end

ease.save_MF3D(); 

