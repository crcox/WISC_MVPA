function ac = ModelContainer(varargin)
% ADLASCONTAINER Expand metavariables to create a structured array, with a
% field for an Model instance.
%
% Input:
% -----
% ModelContainer(expand_struct) : a structure with values that a praticular variable will
% assume over the set of planned experiments.
% ModelContainer(key1,val1,key2,val2,...) : Key value pairs that can be
% coersed into expand_struct.
%
% HYPERBAND special field:
% -----------------------
% A field called HYPERBAND (case sensitive) can be provided, which must
% contain a struct of parameters the correspond to configurations to be
% explored by Hyperband. These variables are not crossed, but are expanded
% with reference to a configID. The configID field will be added to the
% expand struct with the values 1:ncfg.
%
% Output:
% ------
% ac : A structured array, with scalar fields. Each element of the array
% corresponds to a particular experiment configuration.
%
% Chris Cox 07/12/2017
    if nargin == 1 && isstruct(varargin{1})
        expand_struct = varargin{1};
    else
        expand_struct = struct(varargin{:});
    end
    if isfield(expand_struct, 'HYPERPARAMETERS')
        HP = expand_struct.HYPERPARAMETERS;
        if HP.hyperband
            HB = structfun(@pack, rmfield(HP, 'hyperband'), 'UniformOutput', 0);
            HB_fields = fieldnames(HB);
            ncfg = structfun(@numel, HB);
            if all(ncfg == ncfg(1))
                ncfg = ncfg(1);
            elseif all(ncfg == ncfg(1) | ncfg == 1) || all(ncfg == ncfg(end) | ncfg == 1)
                HB_fields_1 = HB_fields(ncfg == 1);
                ncfg = max(ncfg);
                for i = 1:numel(HB_fields_1)
                    f = HB_fields_1{i};
                    HB.(f) = repmat(HB.(f), 1, ncfg);
                end
            else
                error('Hyperband parameter lists must all be of equal length (or length=1).');
            end
            expand_struct = rmfield(expand_struct, 'HYPERPARAMETERS');
            expand_struct.configID = 1:ncfg;
            hbn = fieldnames(HB);
        else
            % If not for use with hyperband, add the parameters back into the
            % main structure. Then drop the HYPERPARAMETERS field.
            HP = structfun(@pack, rmfield(HP, 'hyperband'), 'UniformOutput', 0);
            hpn = fieldnames(HP);
            hbn = {};
            for i = 1:numel(hpn)
                expand_struct.(hpn{i}) = HP.(hpn{i});
            end
            expand_struct = rmfield(expand_struct,'HYPERPARAMETERS');
        end
    end
    fn = fieldnames(expand_struct); n = numel(fn);
    expand_struct_p = structfun(@pack, expand_struct, 'UniformOutput', 0);
    expand_struct_pu = structfun(@safeunique, expand_struct_p, 'UniformOutput', 0);
    for i = 1:n
        if isempty(expand_struct_pu.(fn{i}));
            expand_struct_pu.(fn{i}) = {[]};
        end
    end
    expand_index = struct2cell(structfun(@(x) 1:numel(x), expand_struct_pu, 'UniformOutput', 0));
    B = cell(1,n);
    [B{:}] = ndgrid(expand_index{:});
    for i = 1:n
        if iscell(expand_struct_pu.(fn{i})) % char or cellstr types, or cells containing arrays
            expand_struct.(fn{i}) = expand_struct_pu.(fn{i})(B{i}(:)');
        else % numeric types
            expand_struct.(fn{i}) = num2cell(expand_struct_pu.(fn{i})(B{i}(:)'));
        end
    end
    if isempty(hbn)
        tmp = [struct2cell(expand_struct); {[]}];
        fn = [fn; {'Model'}];
    else
        for i = 1:numel(hbn)
            if iscell(HB.(hbn{i}))
                HB.(hbn{i}) = HB.(hbn{i})(cell2mat(expand_struct.configID));
            else
                HB.(hbn{i}) = num2cell(HB.(hbn{i})(cell2mat(expand_struct.configID)));
            end
        end
        tmp = [struct2cell(expand_struct); struct2cell(HB); {[]}];
        fn = [fn; hbn; {'Model'}];
    end
    expand_cell = [fn(:), tmp(:)]'; % force to be a row vectors
    ac = struct(expand_cell{:});
end
function c = pack(x)
    if ischar(x)
        c = {x};
    else
        c = x;
    end
end
function u = safeunique(x)
    if iscell(x)
        if any(cellfun(@isnumeric, x))
            u = x;
            return
        end
    end
    u = unique(x);
end
            
