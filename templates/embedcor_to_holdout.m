% get predicted coordinates
Cz = embedcor(1).final.Cz;
% get permutation coordinates
Cz_perm = embedcor(1).perm.Cz;
% get cvscheme
cv = meta_tbl.metadata(1).cvind(:,1);
% get target embeddings
dilkina_norms = meta_tbl.metadata(1).targets(2).target;
U = embed_similarity_matrix(dilkina_norms,3);

% initialise output - dimensions x participants x folds
results_all = zeros(length(Cz),3,10);
results_animate = zeros(length(Cz),3,10);
results_inanimate = zeros(length(Cz),3,10);
perms_all = zeros(length(Cz_perm),3,10);
perms_animate = zeros(length(Cz_perm),3,10);
perms_inanimate = zeros(length(Cz_perm),3,10);

% final first
% for each participant
for s = 1:length(Cz)
    % for each dimension
    for d = 1:3
        % for each holdout set
        for ho = 1:10 
            % get target embeddings for that holdout set
            target = U(logical(cv==ho),d);
            % get predicted embeddings for that holdout set
            predicted = Cz{s}(logical(cv==ho),d);
            % correlate
            results_all(s,d,ho) = corr(target,predicted);
            results_animate(s,d,ho) = corr(target(1:5),predicted(1:5));
            results_inanimate(s,d,ho) = corr(target(6:10),predicted(6:10));
        end
    end
end

% then perm
% for each participant/random seed combination
for s = 1:length(Cz_perm)
    % for each dimension
    for d = 1:3
        % for each holdout set
        for ho = 1:10 
            % get target embeddings for that holdout set
            target = U(logical(cv==ho),d);
            % get predicted embeddings for that holdout set
            predicted = Cz_perm{s}(logical(cv==ho),d);
            % correlate
            perms_all(s,d,ho) = corr(target,predicted);
            perms_animate(s,d,ho) = corr(target(1:5),predicted(1:5));
            perms_inanimate(s,d,ho) = corr(target(6:10),predicted(6:10));
        end
    end
end

% average over holdout folds
results_all = mean(results_all,3);
results_animate = mean(results_animate,3);
results_inanimate = mean(results_inanimate,3);
perms_all = mean(perms_all,3);
perms_animate = mean(perms_animate,3);
perms_inanimate = mean(perms_inanimate,3);

% construct permutation distribution.
% perms are currently in the form participant 1, random seed 1; participant
% 2, random seed 1, ... participant 27, random seed 1; participant 1,
% random seed 2... - reshape to get participants x dimensions x random
% seeds
permdist_all = zeros(length(Cz),3,100);
permdist_animate = zeros(length(Cz),3,100);
permdist_inanimate = zeros(length(Cz),3,100);

for seed = 1:100
    permdist_all(:,:,seed) = perms_all(27*(seed-1)+1:27*seed,:);
    permdist_animate(:,:,seed) = perms_animate(27*(seed-1)+1:27*seed,:);
    permdist_inanimate(:,:,seed) = perms_inanimate(27*(seed-1)+1:27*seed,:);
end

% bootstrap by randomly sampling one perm from each participant 10,000
% times
bootstrappedpermdist_all = zeros(length(Cz),3,10000);
bootstrappedpermdist_animate = zeros(length(Cz),3,10000);
bootstrappedpermdist_inanimate = zeros(length(Cz),3,10000);

for i = 1:10000
    % for each participant
    for s = 1:length(Cz)
        % fill in the bootstrapped distribution with random rows of permdist
        bootstrappedpermdist_all(s,:,i) = permdist_all(s,:,randi(100,1));
        bootstrappedpermdist_animate(s,:,i) = permdist_animate(s,:,randi(100,1));
        bootstrappedpermdist_inanimate(s,:,i) = permdist_inanimate(s,:,randi(100,1));        
    end
end

% average over participants to construct the final permutation
% distributions
bootstrappedpermdist_all = squeeze(mean(bootstrappedpermdist_all,1))';
bootstrappedpermdist_animate = squeeze(mean(bootstrappedpermdist_animate,1))';
bootstrappedpermdist_inanimate = squeeze(mean(bootstrappedpermdist_inanimate,1))';

% test whether any permutation distribution differs significantly from 0
bootstrappedpermdist_t_val = zeros(1,9);
bootstrappedpermdist_p_val = zeros(1,9);

[h,p,ci,stats] = ttest(bootstrappedpermdist_all(:,1));
bootstrappedpermdist_p_val(1) = p;
bootstrappedpermdist_t_val(1) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_all(:,2));
bootstrappedpermdist_p_val(2) = p;
bootstrappedpermdist_t_val(2) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_all(:,3));
bootstrappedpermdist_p_val(3) = p;
bootstrappedpermdist_t_val(3) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_animate(:,1));
bootstrappedpermdist_p_val(4) = p;
bootstrappedpermdist_t_val(4) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_animate(:,2));
bootstrappedpermdist_p_val(5) = p;
bootstrappedpermdist_t_val(5) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_animate(:,3));
bootstrappedpermdist_p_val(6) = p;
bootstrappedpermdist_t_val(6) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_inanimate(:,1));
bootstrappedpermdist_p_val(7) = p;
bootstrappedpermdist_t_val(7) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_inanimate(:,2));
bootstrappedpermdist_p_val(8) = p;
bootstrappedpermdist_t_val(8) = stats.tstat;
[h,p,ci,stats] = ttest(bootstrappedpermdist_inanimate(:,3));
bootstrappedpermdist_p_val(9) = p;
bootstrappedpermdist_t_val(9) = stats.tstat;

% some are significantly larger than 0.

% calculate percentile p-values
percentile_p_val = zeros(1,9);
sortedperms = sort(bootstrappedpermdist_all(:,1));
groupmean = mean(results_all(:,1));
b = sum(sortedperms > groupmean);
percentile_p_val(1) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_all(:,2));
groupmean = mean(results_all(:,2));
b = sum(sortedperms > groupmean);
percentile_p_val(2) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_all(:,3));
groupmean = mean(results_all(:,3));
b = sum(sortedperms > groupmean);
percentile_p_val(3) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_animate(:,1));
groupmean = mean(results_animate(:,1));
b = sum(sortedperms > groupmean);
percentile_p_val(4) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_animate(:,2));
groupmean = mean(results_animate(:,2));
b = sum(sortedperms > groupmean);
percentile_p_val(5) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_animate(:,3));
groupmean = mean(results_animate(:,3));
b = sum(sortedperms > groupmean);
percentile_p_val(6) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_inanimate(:,1));
groupmean = mean(results_inanimate(:,1));
b = sum(sortedperms > groupmean);
percentile_p_val(7) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_inanimate(:,2));
groupmean = mean(results_inanimate(:,2));
b = sum(sortedperms > groupmean);
percentile_p_val(8) = (b + 1)/(10000 + 1);
sortedperms = sort(bootstrappedpermdist_inanimate(:,3));
groupmean = mean(results_inanimate(:,3));
b = sum(sortedperms > groupmean);
percentile_p_val(9) = (b + 1)/(10000 + 1);

% plot
figure;
colnames = {'all-D1','all-D2','all-D3','animate-D1','animate-D2','animate-D3','inanimate-D1','inanimate-D2','inanimate-D3'};
b = bar(colnames,mean([results_all,results_animate,results_inanimate]))
ylim([-0.1,1])
stars = {'*','*','*','*','*','','*','','*',};
text(1:9,mean([results_all,results_animate,results_inanimate])+0.01,stars)