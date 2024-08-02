% add useful toolboxes to path 
addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\WISC_MVPA\src')
addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\RSL_correlation_analysis')
addpath('C:\Users\slfri\ownCloud\7T_WISC_MVPA\scripts\RSL_correlation_analysis\m')

%%  Setup (check that you are happy with the whole of this section before running!)

% specify the name of the directory where metadata are found
METADATA_DIR = 'C:\Users\slfri\ownCloud\7T_WISC_MVPA\data';

% specify the name of the directory in which \final and \perm subdirectories
% are found
ANALYSIS_DIR = 'C:\Users\slfri\ownCloud\7T_WISC_MVPA\derivatives\correlation\grOWL\performance';

% specify other parameters. The significance of each is explained in the
% comments of correlation_analysis().
conditions = expand2table(struct( ...
    'TARGET_SOURCE', ["Dilkina_Normalized"], ...
    'AVERAGE_PREDICTED_EMBEDDINGS', [false], ...
    'FULL_RANK_TARGET', [false], ...
    'GrOWL', [true] ...
));
% if you want to run the analysis for multiple targets at once, use this
% syntax, ensuring that the fastest-changing item in the table (i.e. the
% one that flips between configurations the fastest) is first
% conditions = expand2table(struct( ...
%     'TARGET_SOURCE', ["SIFT", "CAFFE7", "CAFFE8", "NEXT6D"], ...
%     'AVERAGE_PREDICTED_EMBEDDINGS', [false], ...
%     'FULL_RANK_TARGET', [false], ...
%     'GrOWL', [true] ...
% ));

% if AVERAGE_PREDICTED_EMBEDDINGS is true, BOOTSTRAP_AVERAGED_EMBEDDINGS
% becomes an option. The default is true when AVERAGE_PREDICTED_EMBEDDINGS
% is true and false when it is false, but it can be switched on and off
% independently too.
% conditions.BOOTSTRAP_AVERAGED_EMBEDDINGS = conditions.AVERAGE_PREDICTED_EMBEDDINGS;
conditions.BOOTSTRAP_AVERAGED_EMBEDDINGS = 1;

% leaving indir options empty means that the script looks for
% final_performance.mat and perm_performance.mat in ANALYSIS_DIR.
indir_options = [];
% Alternatively, this syntax can be used to locate the necessary .mat files
% within file trees (and so to do analyses for mutliple targets at
% once).
% indir_options = [
    % fullfile("SIFT",   "correlation", "performance")
    % fullfile("CAFFE7", "correlation", "performance")
    % fullfile("CAFFE8", "correlation", "performance")
    % fullfile("NEXT6D", "correlation", "performance")
% ];

% setting outdir_options like this means that a results directory will be
% created in ANALYSIS_DIR.
outdir_options = [];
% alternatively, this syntax can be used to save multiple files in
% different locations.
% outdir_options = [
    % fullfile(ANALYSIS_DIR, "SIFT",   "correlation", "low-rank-target", "subject-embeddings")
    % fullfile(ANALYSIS_DIR, "CAFFE7", "correlation", "low-rank-target", "subject-embeddings")
    % fullfile(ANALYSIS_DIR, "CAFFE8", "correlation", "low-rank-target", "subject-embeddings")
    % fullfile(ANALYSIS_DIR, "NEXT6D", "correlation", "low-rank-target", "subject-embeddings")
% ];


%% Open a parallel pool - currently broken!
% Requires parallel computing toolbox
% if isempty(gcp('nocreate'))
%     % pp = parpool(16);
%     pp = parpool(8);
% end

%% run analysis

% start timing the whole script. It can take quite a while!
script_timer = tic();

% onCleanup creates an object that, when destroyed, executes a function -
% in this case, clearing the text progress bar.
cleanupObj = onCleanup(@() clear('textprogressbar'));

% for every target
for i = 1:height(conditions)

    % print a summary of the parameters you have set (as a sanity check)
    c = conditions(i, :);
    disp(table2struct(removevars(c, "BOOTSTRAP_AVERAGED_EMBEDDINGS")));

    %% load data
    % specify directory to load data from
    indir = fullfile(ANALYSIS_DIR,indir_options);
    % print it
    fprintf("Reading data from: %s\n", indir);
    
    % specify directory to write results to
    outdir = fullfile(ANALYSIS_DIR, "results");
    % print it
    fprintf("Writing output to: %s\n", outdir);
    % if the results directory does not exist, make it
    if ~isfolder(outdir)
        mkdir(outdir);
    end
    out_types = ["fullmat", "itemwise", "embedcor"];
    for i = 1:length(out_types)
        if ~isfolder(fullfile(outdir, out_types(i)));
            mkdir(fullfile(outdir, out_types(i)));
        end
    end

    % start timing how long it takes to load the data and format it
    tic;
    fprintf('Loading data ');

    % load final results into a variable called final
    final = getfield(load(fullfile(indir, 'final_performance.mat')), 'Tallcv');
    % load perm results into a variable called perm
    perm = getfield(load(fullfile(indir, 'perm_performance.mat')), 'Tallcv');

    %% Convert cell array of character vectors into strings
    final.LambdaSeq = repmat("inf", height(final), 1);
    perm.LambdaSeq = repmat("inf", height(perm), 1);

    % list variables to be converted into strings
    vars = ["filters", "FiltersToApplyBeforeEmbedding", "target_label", ...
        "target_type", "sim_source", "sim_metric", "data", "data_varname", ...
        "metadata", "metadata_varname", "normalize_wrt", "normalize_data", ...
        "normalize_target", "LambdaSeq"];

        % convert those variables into strings
    for v = vars
        final.(v) = cellstr2string(final.(v));
        perm.(v) = cellstr2string(perm.(v));
    end

    % bug fix - if for any reason you have set the analysis to output
    % results.mat to anywhere other than the directory from which the
    % analysis was run (e.g. /staging), there will be no params.json handy
    % which affects the way that load_from_condor works. Including these
    % lines will make this script run unimpeded.
    % final.filters = repmat({["rowfilter","colfilter"]},size(final,1),1);
    % final.FiltersToApplyBeforeEmbedding = repmat(cell(1,1),size(final,1),1);
    % perm.filters = repmat({["rowfilter","colfilter"]},size(perm,1),1);
    % perm.FiltersToApplyBeforeEmbedding = repmat(cell(1,1),size(perm,1),1);
    

    %% load metadata

    % load
    tmp = load([METADATA_DIR,'\metadata.mat'], 'metadata');
    % create a table containing one copy of the metadata per participant
    % plus the name of the variable containing the metadata
    meta_tbl = table(tmp.metadata, "metadata", 'VariableNames', ["metadata", "metadata_varname"]);
    
    % finish timing data loading and formatting
    fprintf('(%.2f s)\n', toc);


    %% Main analysis
    [fullmat, itemwise, embedcor] = ...
      correlation_analysis( ...
        final, perm, meta_tbl, ...
        'window_type','None', ...
        'category_labels', ["animate","inanimate"], ...
        'full_rank_target', c.FULL_RANK_TARGET, ...
        'average_predicted_embeddings', c.AVERAGE_PREDICTED_EMBEDDINGS, ...
        'bootstrap_averaged_embeddings', c.BOOTSTRAP_AVERAGED_EMBEDDINGS, ...
        'cvscheme', 1, ...
        'scale_singular_vectors', true);
        

    %% Write text
    tic;
    textprogressbar(sprintf('%36s', 'Writing correlations to text: '));
    textprogressbar(0);
    writetable(fullmat(end).final, fullfile(outdir, "fullmat", "final.csv"));
    textprogressbar((1/6)*100);
    writetable(fullmat(end).perm, fullfile(outdir, "fullmat", "perm.csv"));
    textprogressbar((2/6)*100);
    writetable(itemwise(end).final, fullfile(outdir, "itemwise", "final.csv"));
    textprogressbar((3/6)*100);
    writetable(itemwise(end).perm, fullfile(outdir, "itemwise", "perm.csv"));
    textprogressbar((4/6)*100);
    writetable(embedcor(end).final, fullfile(outdir, "embedcor", "final.csv"));
    textprogressbar((5/6)*100);
    writetable(embedcor(end).perm, fullfile(outdir, "embedcor", "perm.csv"));
    textprogressbar((6/6)*100);
    textprogressbar(sprintf(' done (%.2f s)', toc));


    %% Save mat
    tic;
    textprogressbar(sprintf('%36s', 'Saving correlations to .mat: '));
    save(fullfile(outdir, "fullmat.mat"), 'fullmat');
    textprogressbar((1/3)*100);
    save(fullfile(outdir, 'itemwise.mat'), 'itemwise');
    textprogressbar((2/3)*100);
    save(fullfile(outdir, 'embedcor.mat'), 'embedcor');
    textprogressbar((3/3)*100);
    textprogressbar(sprintf(' done (%.2f s)', toc));

    fprintf("\n");
    fprintf("----");

end

toc(script_timer)

%% Close parallel pool
% delete(gcp('nocreate'));
