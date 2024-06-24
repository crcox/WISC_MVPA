function load_tuning_performance(path)
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\src')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\util')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\dependencies\jsonlab')
    addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts') 
    
    Tallcv = load_from_condor(path,'skip_large_matrices',true); 
    save([fileparts(path),'\tune_performance.mat'],'Tallcv');
    writetable(Tallcv,[fileparts(path),'\tune_performance.csv']);
end
