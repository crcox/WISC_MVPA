classdef Subject
    %SUBJECT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        filename
        subject
        label
        data
        cvscheme
        rowfilters
        colfilters
        targets
        coords
        permutations
    end
    
    methods
        % This function should be more stringent. For example, check that
        % required fields are set and check that they are set to sane
        % values.
        function obj = Subject(varargin)
            if nargin == 0
                return
            elseif nargin == 1
                s = varargin{1};
            else
                s = struct(varargin{:});
            end
            fn = fieldnames(s);
            for i = 1:numel(fn)
                disp(fn{i})
                obj.(fn{i}) = s.(fn{i});
            end
        end
        
        function rf = getRowFilter(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});
            
            if ~isempty(p.Results.include)
                RF = selectbyfield(obj.rowfilters, 'label', p.Results.include, 'dimension', 1);
            elseif ~isempty(p.Results.exclude)
                labels = {obj.rowfilters.label};
                z = ~ismember(labels, p.Results.exclude);
                RF = selectbyfield(obj.rowfilters, 'label', labels(z), 'dimension', 1);
            else
                RF = obj.rowfilters;
            end
            if numel(RF) == 0
                RF(1).filter = true(1, size(obj.data, 1));
            end
            rf = all(cat(1, RF.filter),1);
        end
        
        function cf = getColFilter(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});
            
            if ~isempty(p.Results.include)
                CF = selectbyfield(obj.colfilters, 'label', p.Results.include, 'dimension', 2);
            elseif ~isempty(p.Results.exclude)
                labels = {obj.colfilters.label};
                z = ~ismember(labels, p.Results.exclude);
                CF = selectbyfield(obj.colfilters, 'label', labels(z), 'dimension', 2);
            else
                CF = obj.colfilters;
            end
            if numel(CF) == 0
                CF(1).filter = true(1, size(obj.data, 2));
            end
            cf = all(cat(1, CF.filter),1);
        end
        
        function [rf,cf] = getFilters(obj, varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});
            
            rf = obj.getRowFilter(varargin{:});
            cf = obj.getColFilter(varargin{:});
        end
        
        function c = getCoords(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addOptional(p, 'orientation', '', @ischar);
            addOptional(p, 'coordsystem', '', @ischar);
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            parse(p, obj, varargin{:});
            
            if p.Results.unfiltered
                cf = true(1, size(obj.data,2));
            else
                if ~isempty(p.Results.include);
                    cf = obj.getColFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    cf = obj.getColFilter('exclude',p.Results.exclude);
                else
                    cf = obj.getColFilter();
                end
            end
            if ~isempty(p.Results.orientation)
                x = selectbyfield(obj.coords,'orientation', p.Results.orientation);
            else
                x = obj.coords;
            end
            if ~isempty(p.Results.coordsystem)
                x.(p.Results.coordsystem) = x.(p.Results.coordsystem)(cf, :);
            else
                for i = 1:numel(x)
                    if isfield(x(i), 'ind') && ~isempty(x(i).ind), x(i).ind = x(i).ind(cf,:); end
                    if isfield(x(i), 'ijk') && ~isempty(x(i).ijk), x(i).ijk = x(i).ijk(cf,:); end
                    if isfield(x(i), 'xyz') && ~isempty(x(i).xyz), x(i).xyz = x(i).xyz(cf,:); end
                end
            end
            if p.Results.simplify && ~isempty(p.Results.coordsystem) && (numel(x) == 1)
                c = x(1).(p.Results.coordsystem);
            else
                c = x;
            end
        end
        
        function x = getCVScheme(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});
            
            if p.Results.unfiltered
                x = obj.cvscheme;
            else
                if ~isempty(p.Results.include);
                    rf = obj.getRowFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    rf = obj.getRowFilter('exclude',p.Results.exclude);
                else
                    rf = obj.getRowFilter();
                end
                x = obj.cvscheme(rf);
            end
        end        
        function x = getData(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            parse(p, obj, varargin{:});
            
            if p.Results.unfiltered
                x = obj.data;
            else
                if ~isempty(p.Results.include);
                    [rf, cf] = obj.getFilters('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    [rf, cf] = obj.getFilters('exclude',p.Results.exclude);
                else
                    [rf, cf] = obj.getFilters();
                end
                x = obj.data(rf, cf);
            end
        end
        
        function c = getTargets(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(isnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(isnumeric(x) == [2,1,0]));
            parse(p, obj, varargin{:});
            
            if p.Results.unfiltered
                rf = true(1, size(obj.data,1));
            else
                if ~isempty(p.Results.include);
                    rf = obj.getRowFilter('include',p.Results.include);
                elseif ~isempty(p.Results.exclude);
                    rf = obj.getRowFilter('exclude',p.Results.exclude);
                else
                    rf = obj.getRowFilter();
                end
            end

            x = obj.targets;
            for i = 1:numel(x)
                if iscell(x(i).type)
                    for j = 1:numel(x(i).type)
                        if strcmpi(x(i).type{j},'similarity')
                            x(i).target{j} = x(i).target{j}(rf,rf);
                        else
                            x(i).target{j} = x(i).target{j}(rf,:);
                        end
                    end
                else
                    if strcmpi(x(i).type,'similarity')
                        x(i).target = x(i).target(rf,rf);
                    else
                        x(i).target = x(i).target(rf,:);
                    end
                end
            end
            if p.Results.simplify && (numel(x) == 1)
                if iscell(x(1).target)
                    % This is poorly implemented, but simplify can take on
                    % either a logical value or a numeric value of 0, 1,
                    % *or 2*. If 2, then we select from "the history". By
                    % this, I am refering to the fact that similarity
                    % matrices tagged with type 'similarity' might be
                    % converted to embeddings. In these cases, multiple
                    % target values will be stored in a cell array, with
                    % most recent version (the embedding) being the first
                    % element and the prior version (the similarity matrix)
                    % being the second.
                    % The overloading of the simplicity parameter is
                    % further abused by exploiting that true behaves the
                    % same way as 1 in this context. This means if simplify
                    % == true or simplify == 1, the more recent version of
                    % the target (the embedding) will be selected, and if
                    % simplify == 2 the similarity matrix will be selected.
                    c = x(1).target{p.Results.simplify};
                else
                    c = x(1).target;
                end
            else
                c = x;
            end
        end
        
        function pix = getPermutationIndex(obj, RandomSeed)
            z = obj.permutations.RandomSeed == RandomSeed;
            pix = obj.permutations.index(:,z);
        end
        
        function tp = getPermutedTargets(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addParameter(p, 'unfiltered', false, @(x) islogical(x) || any(isnumeric(x) == [1,0]));
            addParameter(p, 'include', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'exclude', {}, @(x) iscell(x) || ischar(x) );
            addParameter(p, 'simplify', false, @(x) islogical(x) || any(isnumeric(x) == [2,1,0]));
            addParameter(p, 'RandomSeed', 0, @(x) isnumeric(x) || isscalar(x));
            parse(p, obj, varargin{:});
            
            t = obj.getTargets(...
                'unfiltered', true, ...
                'include', p.Results.include, ...
                'exclude', p.Results.exclude, ...
                'simplify', p.Results.simplify);
            pix = obj.getPermutationIndex(p.Results.RandomSeed);
            if p.Results.unfiltered
                tp = t(pix);
            else
                rf = obj.getRowFilter();
                rfp = rf(pix);
                tp = t(pix);
                tp = tp(rfp);
            end
        end
        
        function obj = generateEmbeddings(obj, tau, varargin)
            p = inputParser();
            addRequired(p, 'obj', @(x) isa(x, 'Subject'));
            addRequired(p, 'tau', @(x) isnumeric(x) && isscalar(x));
            addParameter(p, 'ExtendEmbedding', false, @(x) islogical(x) || any(asnumeric(x) == [1,0]));
            parse(p, obj, tau, varargin{:});
            for i = 1:numel(obj.targets)
                if strcmpi(obj.targets(i).type, 'similarity');
                    S = obj.getTargets(i,'unfiltered',true,'simplify',true);
                    [C, r] = sqrt_truncate_r(S, tau);
                    fprintf('S decomposed into %d dimensions (tau=%.2f)\n', r, tau);
                    
                    nexamples = size(obj.getData('unfiltered', true), 1);
                    if size(C,1) < nexamples
                        remainder = rem(nexamples, size(C,1));
                        if remainder > 0
                            error('C has fewer rows than X, and number in X is not evenly divisible by number in C.');
                        else
                            repeatntimes = size(nexamples,1) ./ size(C,1);
                            warning('C has fewer rows than X, and number in X is evenly divisible by number in C. Repeating C %d times to match.', repeatntimes);
                            C = repmat(C, repeatntimes, 1);
                        end
                    end
                    obj.targets(i).target = {C,S};
                    obj.targets(i).type = {'embedding','similarity'};
                end
            end
        end
    end
end

