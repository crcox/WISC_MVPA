function Tallcv = load_from_condor(rootdir, varargin)
    p = inputParser;
    addRequired(p, 'rootdir');
    addParameter(p, 'skip_large_matrices', false, @islogical);
    parse(p, rootdir, varargin{:});
    if p.Results.skip_large_matrices
        SKIP = {'Uz','Cz','nz_rows','coords'};
    else
        SKIP = {};
    end
    x = dir(rootdir);
    z = [x.isdir] & cellfun(@(x) ~isempty(regexp(x, '^\d+$', 'start')), {x.name});
    jobdirs = x(z);
    njobs = length(jobdirs);
    fmt = sprintf('%%0%dd', strlength(jobdirs(1).name));
    for i = 1:njobs
        jobdir = sprintf(fmt, i - 1);
        if exist(fullfile(rootdir,jobdir,'results.mat'),'file')
            tmp = load(fullfile(rootdir,jobdir,'results.mat'),'results');
            if ~isempty(SKIP)
                tmp.results = rmfield(tmp.results, SKIP);
            end
            [tmp.results.jobdir] = deal(i - 1);
            if exist(fullfile(rootdir,jobdir,'params.json'),'file')
                j = loadjson(fullfile(rootdir,jobdir,'params.json'));
                [tmp.results.LambdaSeq] = deal(string(j.LambdaSeq));
                [tmp.results.filters] = deal(j.filters);
                if isfield(j, 'FiltersToApplyBeforeEmbedding')
                    [tmp.results.FiltersToApplyBeforeEmbedding] = deal(j.FiltersToApplyBeforeEmbedding);
                else
                    [tmp.results.FiltersToApplyBeforeEmbedding] = deal({});
                end
            end
            if i == 1
                R = structfun(@(x) [], tmp.results(1), 'UniformOutput', 0);
                R = repmat(R,numel(tmp.results),njobs,1);
            end
            R(:,i) = tmp.results(:);
        end
    end
    z = ~cellfun(@isempty, {R.jobdir});
    results = R(z);
    clear tmp R

    % Compose a table of the data
    Tallcv = struct2table(results(:));
end
