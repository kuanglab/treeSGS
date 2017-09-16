clc, clear;
close all;
format compact;

addpath('../exonTLFL/');
addpath(genpath('../gflasso/SPG_Multi_Graph'));
%% load and process data
load Hapmap_CNV_processed2;

X = data;
samples = ThePop;
[m, n] = size(X);

%% PCA plot
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

save allPop_info mx group poplist pop_num pop_ix;
%% parameters
figure(1);
figure(2);
iii = 1;
for k = 120
for ratio = [0.5 0.6 0.9]
for pcut = [0.00005 0.0001]
for sigma = [20 40]

gamma = 0;
lambda = 1;
Vnorm = 2;

maxiter = 100;
[m, n] = size(X);

%%
load(['TBL6/TBL6_result/TreeLassoothertest6_allpop_k' num2str(k) '_ratio' num2str(ratio) ...
    '_pcut' num2str(pcut) '_sigma' num2str(sigma) ...
    '_lambda' num2str(lambda) '_gamma' num2str(gamma) ...
    '_Vnorm' num2str(Vnorm) '.mat']);

% re-arrange
mx_01 = zeros(k, pop_num);
for i = 1:pop_num
    mx_01(:,  i) = sum(V_sps(:, pop_ix{i}), 2) > 0;
end;

newix = 1:k;
newmx_01 = mx_01;
for i = pop_num:-1:1
    [~, IX] = sort(newmx_01(:,i), 'descend');
    newmx_01 = newmx_01(IX, :);
    newix = newix(IX);
end;

V_sps_new = V_sps(newix,:);
U_sps_new = U_sps(:,newix);


% figure;
% for i = 1:k
%     subplot (5,k/5,i);
%     plot(U_sps_new(:,i));
%     title(i);
% end;

figure(1);
subplot(2,2,iii);
imagesc(V_sps_new);
set(gca,'XTick', cellfun(@mean, pop_ix), 'xaxisLocation','top');
set(gca,'XTickLabel', poplist, 'xaxisLocation','top');
hold on;
for i = 1:pop_num
    line([pop_ix{i}(end)+0.5, pop_ix{i}(end)+0.5], [0.5, k+0.5], 'Color', 'r');
end;
hold off;
colorbar();
title(['k:' num2str(k) ', ratio:' num2str(ratio) ',sigma:' num2str(sigma) ' ,p:' num2str(pcut)]);


figure(2);
subplot(2,4,iii*2-1);
plot(Cov_sps);
title(['cov, sigma:' num2str(sigma) ' ,p:' num2str(pcut)]);
subplot(2,4,iii*2);
plot(Obj_sps(4:end,1));
title(['obj, sigma:' num2str(sigma) ' ,p:' num2str(pcut)]);

iii = iii + 1;
end;
end;
end;
end;