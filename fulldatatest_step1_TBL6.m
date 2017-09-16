clc, clear;
close all;
format compact;

addpath('../exonTLFL/');
addpath(genpath('../gflasso/SPG_Multi_Graph'));
addpath(genpath('../SLEP_package_4.1/SLEP'));

%% load and process data
load Hapmap_CNV_processed2;

X = data;
samples = ThePop;
[m, n] = size(X);

%% PCA
% [~, score, latent] = pca(X');
% ft = 10;
% imagesc(score(:,1:ft));
% set(gca, 'YAxisLocation', 'right');
% hold on;
% line([0, ft+0.5], [511.5, 511.5], 'Color', 'r');
% line([0, ft+0.5], [588.5, 588.5], 'Color', 'r');
% line([0, ft+0.5], [842.5, 842.5], 'Color', 'r');
% line([0, ft+0.5], [1093.5, 1093.5], 'Color', 'r');
% hold off;
% title('PCA of X');

%% phylogenetic info
% pop seqence in data
poplist = {'ASW', 'LWK', 'MKK', 'YRI', 'MXL', 'CHB', 'CHD', 'JPT', 'CEU', 'TSI', 'GIH'};
pop_num = length(poplist);

pop_ix = cell(pop_num,1);
for i = 1:pop_num
    pop_ix{i} = find(strcmp(ThePop(:,3), poplist{i}));
end;

% group
group = zeros(pop_num, n);
for i = 1:pop_num
    group(i, pop_ix{i}) = 1;
end;

% graph
mx = zeros(pop_num*2-1);
mx([2 3], 12) = 1;
mx([4 12], 13) = 1;
mx([1 13], 14) = 1;
mx([6 7], 15) = 1;
mx([8 15], 16) = 1;
mx([9 10], 17) = 1;
mx([5 17], 18) = 1;
mx([11 18], 19) = 1;
mx([16 19], 20) = 1;
mx([14 20], 21) = 1;

%% parameters
for k = [80:10:120]
for ratio = [0.8]
for pcut = [0.0001 0.0005 0.001 0.005]
for sigma = [20]
disp([k ratio pcut sigma]);
lambda = 1;
gamma = 0;
Vnorm = 2;

maxiter = 100;
[m, n] = size(X);

C = zeros(m-1, m);

%% tree based
tStart = tic;


[U_sps, V_sps, Cov_sps, Obj_sps, U_scale_his_sps] = ...
    FL_exon_TreebaseLasso6(X, gamma, lambda, k, C, maxiter, [], [], group, mx, Vnorm, ratio, pcut, sigma);

tElapsed = toc(tStart);
numIter = size(Cov_sps,1);

save(['TBL6/TBL6_result/TreeLassoothertest6_allpop_k' num2str(k) '_ratio' num2str(ratio) ...
    '_pcut' num2str(pcut) '_sigma' num2str(sigma) ...
    '_lambda' num2str(lambda) '_gamma' num2str(gamma) ...
    '_Vnorm' num2str(Vnorm) '_runtime.mat'],...
    'numIter', 'tElapsed');


% save(['TBL6/TBL6_result/TreeLassoothertest6_allpop_k' num2str(k) '_ratio' num2str(ratio) ...
%     '_pcut' num2str(pcut) '_sigma' num2str(sigma) ...
%     '_lambda' num2str(lambda) '_gamma' num2str(gamma) ...
%     '_Vnorm' num2str(Vnorm) '.mat'],...
%     'U_sps', 'V_sps', 'Cov_sps', 'Obj_sps', 'U_scale_his_sps');
end;
end;
end;
end;