classdef Searchlight
    %SEARCHLIGHT Summary of this class goes here
    %   Detailed explanation goes here

    properties
        radius
        nvox
        SL
        error_map1
        error_map2
        trainingFilter
    end

    properties ( Access = public, Hidden = true )
        s
        max_iter = 100000
        tol = 1e-8
        ModelHasBiasUnit = false
        glmnetParallel = false
        glmnetCV
        currentSearchlight = 0
        failedSearchlights
        EMPTY = 1
    end

    methods
        function obj = Searchlight(radius,coords,cvind,trainingFilter,ModelHasBiasUnit,opts)
            if (nargin == 0)
                return
            end
            if (nargin < 4), opts = struct(); end
            obj.radius = radius;
            obj.SL = obj.generate_searchlights(coords, radius);
            obj.ModelHasBiasUnit = ModelHasBiasUnit;

            obj.nvox = size(coords, 1);
            obj.trainingFilter = trainingFilter;

            [~,~,cvind_adj] = unique(cvind(trainingFilter));
            obj.glmnetCV = cvind_adj;

            obj.error_map1 = nan(obj.nvox, 1);
            obj.error_map2 = nan(obj.nvox, 1);
            obj.failedSearchlights = false(obj.nvox, 1);

            fn = fieldnames(opts);
            for i = 1:numel(fn)
                obj.(fn{i}) = opts.(fn{i});
            end
            obj.setRandomSeed(0);
            obj.EMPTY = 0;
        end

        function obj = setRandomSeed(obj, seed)
            obj.s = RandStream('mt19937ar','Seed',seed);
        end

        function SL = generate_searchlights(~,xyz,radius)
            D = squareform(pdist(xyz));
            SL = cell(size(D,1),1);
            for j = 1:size(D,1)
                SL{j} = find(D(j,:) <= radius);
            end
        end

        function obj = computeInformationMap(obj,X,C)
            [Xtrain,Ctrain] = obj.getTrainingData(X,C,'TransposeX',false);
            [Xtest,Ctest] = obj.getTestingData(X,C,'TransposeX',false);
            [done,failed] = deal(0);
            for i = 1:numel(obj.SL)
                nchar = fprintf('Progress: %d done, %d failed, out of %d', done, failed, numel(obj.SL));
                obj.currentSearchlight = i;
                g = obj.SL{i};
                [glmnetOpts, FAILED] = obj.tune(Xtrain(:,g),Ctrain);
                if FAILED
                    obj.failedSearchlights(i) = true;
                    failed = failed + 1;
                    eraser = repmat('\b',1,nchar);
                    fprintf(eraser);
                    continue
                end
                glmnetModel = obj.train(Xtrain(:,g),Ctrain,glmnetOpts);
                obj.error_map2(i) = obj.test(Xtrain(:,g),Ctrain,glmnetModel);
                obj.error_map1(i) = obj.test(Xtest(:,g),Ctest,glmnetModel);
                done = done + 1;
                eraser = repmat('\b',1,nchar);
                fprintf(eraser);
            end
        end

        function [opts_tuned, return_code] = tune(obj,Xt,Ct)
            cvn = max(obj.glmnetCV);
            try
                opts_cv = glmnetSet( ...
                    struct( ...
                        'mtype','grouped', ...
                        'alpha',1, ...
                        'intr',obj.ModelHasBiasUnit));
                fitobj_cv = cvglmnet(Xt, Ct, 'mgaussian', opts_cv, 'mse', cvn, obj.glmnetCV, obj.glmnetParallel);
                opts_tuned = glmnetSet( ...
                    struct( ...
                        'mtype','grouped', ...
                        'alpha',1, ...
                        'intr',obj.ModelHasBiasUnit, ...
                        'lambda',fitobj_cv.lambda_min));
                return_code = 0;
            catch
                opts_tuned = [];
                return_code = 1;
            end
        end

        function glmnetModel = train(~,Xt,Ct,opts)
            glmnetModel = glmnet(Xt, Ct, 'mgaussian', opts);
        end

        function loss = test(~,X, C, glmnetModel)
            Cz = glmnetPredict(glmnetModel, X);
            loss = nrsa_loss(C,Cz);
        end

        function [X,C] = getTrainingData(obj,Xa,Ca,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            if p.Results.transposeX
                X = Xa(obj.trainingFilter,:)';
            else
                X = Xa(obj.trainingFilter,:);
            end
            C = Ca(obj.trainingFilter,:);
        end

        function [X,C] = getTestingData(obj,Xa,Ca,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addParameter(p, 'transposeX', false, @(x) islogical(x) || x==1 || x==0);
            parse(p, obj, varargin{:});

            if p.Results.transposeX
                X = Xa(~obj.trainingFilter,:)';
            else
                X = Xa(~obj.trainingFilter,:);
            end
            C = Ca(~obj.trainingFilter,:);
        end
        
        function printresults(obj,varargin)
            p = inputParser();
            addRequired(p, 'obj');
            addOptional(p, 'outputcontrol', 'bycv', @ischar);
            addOptional(p, 'bysubject', 0, @isnumeric);
            addOptional(p, 'subjects', 1, @isnumeric);
            parse(p, obj, varargin{:});

            cvix = p.Results.cvindex;
            subj = p.Results.subjects;
            switch lower(p.Results.outputcontrol)
                case 'header'
                    fprintf('%8s%8s%8s%8s%8s\n', 'subj','cvindex','radius','failures','nvox');
                case 'bysubject'
                    fprintf('%8d%8d%8d%8d%8d\n', ...
                        subj,cvix,obj.radius,nnz(obj.failedSearchlights),obj.nvox);

            end

        end

    end
end

