%% 
ease.get_MF3D(); 
neuron.orderROIs('snr'); 
neuron0 = neuron.copy(); 
RSS0 = neuron.compute_rss(); 
A0 = neuron.A; 
C0 = neuron.C; 
Craw0 = neuron.C_raw; 

K = size(neuron.A, 2); 
maxIter = 20; 

%% run multiple iterations
% with weighting 
with_em = true; 
RSS_withEM = zeros(1, maxIter); 
ccA_withEM = zeros(K, maxIter);
ccC_withEM = zeros(K, maxIter); 
ccCraw_withEM = zeros(K, maxIter); 
neuron = neuron0.copy(); 
for m=1:maxIter 
    neuron.hals([]); 
    RSS_withEM(m) = neuron.compute_rss(); 
    ccA_withEM(:, m) = sum(A0.*neuron.A, 1) ./ sqrt(sum(A0.^2,1).*sum(neuron.A.^2, 1)); 
    ccC_withEM(:,m) = sum(C0.*neuron.C, 2) ./ sqrt(sum(C0.^2,2).*sum(neuron.C.^2, 2));     
    ccCraw_withEM(:,m) = sum(Craw0.*neuron.C_raw, 2) ./ sqrt(sum(Craw0.^2,2).*sum(neuron.C_raw.^2, 2));     
    disp(m); 
end 
neuron_withEM = neuron.copy(); 
neuron = neuron0.copy(); 

% without EM info 
with_em = false; 
RSS_noEM = zeros(1, maxIter); 
ccA_noEM = zeros(K, maxIter);
ccC_noEM = zeros(K, maxIter); 
for m=1:maxIter 
    neuron.hals([], false, with_em); 
    RSS_noEM(m) = neuron.compute_rss(); 
    ccA_noEM(:, m) = sum(A0.*neuron.A, 1) ./ sqrt(sum(A0.^2,1).*sum(neuron.A.^2, 1)); 
    ccC_noEM(:,m) = sum(C0.*neuron.C, 2) ./ sqrt(sum(C0.^2,2).*sum(neuron.C.^2, 2));     
    ccCraw_noEM(:,m) = sum(Craw0.*neuron.C_raw, 2) ./ sqrt(sum(Craw0.^2,2).*sum(neuron.C_raw.^2, 2));     
    disp(m); 
end 
neuron_noEM = neuron.copy(); 
neuron = neuron0.copy(); 

% without weighting 
weight_em = false; 
RSS_noweight = zeros(1, maxIter); 
ccA_noweight = zeros(K, maxIter);
ccC_noweight = zeros(K, maxIter); 
for m=1:maxIter 
    neuron.hals([], weight_em); 
    RSS_noweight(m) = neuron.compute_rss(); 
    ccA_noweight(:, m) = sum(A0.*neuron.A, 1) ./ sqrt(sum(A0.^2,1).*sum(neuron.A.^2, 1)); 
    ccC_noweight(:,m) = sum(C0.*neuron.C, 2) ./ sqrt(sum(C0.^2,2).*sum(neuron.C.^2, 2));     
    ccCraw_noweight(:,m) = sum(Craw0.*neuron.C_raw, 2) ./ sqrt(sum(Craw0.^2,2).*sum(neuron.C_raw.^2, 2));     
    disp(m); 
end 
neuron_noweight = neuron.copy(); 
neuron = neuron0.copy(); 

save EMinfo_results RSS* ccA* ccC* ccCraw* neuron_noEM neuron_withEM neuron0 neuron_noweight

%% find an example neuron 
example_em_id = 73430227; 
example_id = find(neuron.get_em_ids(1:K, EM_info) == example_em_id); 

neuron_noEM.showNeuron(example_id); 
neuron_noweight.showNeuron(example_id); 
neuron_withEM.showNeuron(example_id); 
neuron.showNeuron(example_id); 

%% 
figure; 
plot(RSS_withEM/RSS0, '-o'); 
hold on; 
plot(RSS_noweight/RSS0, '-v'); 
plot(RSS_noEM/RSS0, '-s'); 
legend('with support & weighting', 'with support, no weighting', ...
'no support, no weighting', 'orientation', 'vertical'); 

%% spatial A
vlim = [0, 1]; 
figure; 
subplot(131); 
imagesc(ccA_withEM, vlim); 
hold on; 
axis off; 
title('EASE'); 
subplot(132); 
imagesc(ccA_noweight, vlim);
axis off; 
title('no weighting'); 

subplot(133); 
imagesc(ccA_noEM, vlim); 
pos = get(gca, 'position'); 
colorbar; 
set(gca, 'position', pos); 
axis off; 
title('no EM'); 

%% temporal, C 
vlim = [0, 1]; 
figure; 
subplot(131); 
imagesc(ccC_withEM, vlim); 
hold on; 
axis off; 
title('EASE'); 
subplot(132); 
imagesc(ccC_noweight, vlim);
axis off; 
title('no weighting'); 

subplot(133); 
imagesc(ccC_noEM, vlim); 
pos = get(gca, 'position'); 
colorbar; 
set(gca, 'position', pos); 
axis off; 
title('no EM'); 

%% temporal, Craw 
vlim = [0, 1]; 
figure; 
subplot(131); 
imagesc(ccCraw_withEM, vlim); 
hold on; 
axis off; 
title('EASE'); 
subplot(132); 
imagesc(ccCraw_noweight, vlim);
axis off; 
title('no weighting'); 

subplot(133); 
imagesc(ccCraw_noEM, vlim); 
pos = get(gca, 'position'); 
colorbar; 
set(gca, 'position', pos); 
axis off; 
title('no EM'); 


%% 
figure; 
plot(ccC_withEM(:,20), '-b'); 
hold on; 
plot(ccC_noweight(:,20), '-r'); 
plot(ccC_noEM(:,20), '-g'); 
legend('with support & weighting', 'with support, no weighting', ...
'no support, no weighting', 'location', 'southwest'); 

xlabel('cell #');
ylabel('similarity(ci, ci0)')
ylim([0, 1.01]); 










