%% save all neurons 
em_ids = cell2mat(neuron.match_status.em_ids);
tmp_folder = fullfile(output_figs, 'neurons');
if ~exist(tmp_folder, 'dir')
    mkdir(tmp_folder);
else
    delete(tmp_folder);
    mkdir(tmp_folder);
end
for m=1:length(em_ids)
    neuron.showNeuron(m);
    em_id = em_ids(m);
    export_fig(gcf, fullfile(output_figs, 'neurons', sprintf('neuron_%d_%d.pdf', m, em_id)));
    export_fig(gcf, fullfile(output_figs, 'neurons', sprintf('neuron_%d_%d.fig', m, em_id)));
    pause(0.1);
    close;
end