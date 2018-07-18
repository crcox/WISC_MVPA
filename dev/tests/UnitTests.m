classdef UnitTests < matlab.unittest.TestCase
    properties
        OMIT
        metadata
        params
        Y
        X
    end
    methods (TestClassSetup)
        function LoadData(testCase)
            addpath('src');
            testCase.OMIT = 1;
            load('testdata/metadata.mat');
            testCase.metadata = metadata; %#ok<CPROP>
            params = loadjson('testdata/000/params.json'); %#ok<PROP>
            params = init_opts(params); %#ok<PROP>
            testCase.params = params; %#ok<PROP>
            load('testdata/alldata.mat', 'X');
            testCase.X = X; %#ok<CPROP>
%             Y = {testCase.metadata.TrueFaces}';
        end
    end

    %% Test Method Block
    methods (Test)
        % Includes unit test functions
        function testCommonSpaceDims(testCase) 
          [x,y,z] = ndgrid(1:10,1:10,1:10);
          xyz = [x(:), y(:), z(:)];
          [I,J,K] = defineCommonIndexSpace(xyz,1,[6,5,4],[3,2,2]);
          testCase.assertEqual([I,J,K], [12,11,12]);
				end
				
        function testNonsmoothEval00(testCase)
				% All group 1
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 1; 
					W = [1,1,0,0]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sqrt(1+1))
				end
				
        function testNonsmoothEval01(testCase)
				% All group 2
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 1; 
					W = [0,0,1,1]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sqrt(1+1))
				end
					
        function testNonsmoothEval02(testCase)
				% Across groups
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 1; 
					W = [1,0,0,1]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sqrt(1)+sqrt(1))
				end
				
        function testNonsmoothEval03(testCase)
				% 2 subjects, group 1 only
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 1; 
					W = [1,1,0,0]';
					W = [W,W];
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sqrt(1+1+1+1))
				end
        function testNonsmoothEval10(testCase)
				% All group 1
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 0; 
					W = [1,1,0,0]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sum(abs(W(:))))
				end
				
        function testNonsmoothEval11(testCase)
				% All group 2
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 0; 
					W = [0,0,1,1]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sum(abs(W(:))))
				end
					
        function testNonsmoothEval12(testCase)
				% Across groups
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 0; 
					W = [1,0,0,1]';
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sum(abs(W(:))))
				end
				
        function testNonsmoothEval13(testCase)
				% 2 subjects, group 1 only
					lambda = 1;
					group_arr = [1,2;3,4];
					alpha = 0; 
					W = [1,1,0,0]';
					W = [W,W];
					Fz = nonsmooth_eval(W, lambda, alpha, group_arr);
					testCase.assertEqual(Fz, sum(abs(W(:))))
				end
				function testSOSLassoShrinkLogistic_0(testCase)
					W = [1,1,1,0,0; 1,1,1,1,0; 1,1,1,1,1; 0,0,0,0,0];
					[p,~] = size(W);
					alpha = 1;
					lam = 2;
					G = (1:p)';
					groups = (1:p);
					y = soslasso_shrink_logistic(W,G,groups,lam,alpha);
					nonzero = sqrt(sum(W.^2,2)) > lam;
					y_nonzero = sqrt(sum(y.^2,2)) > 0;
					testCase.assertEqual(nonzero,y_nonzero);
				end
				function testSOSLassoShrinkLogistic_1(testCase)
					W = [1,2,0,0,0; 1,0,2,0,0; 1,0,0,2,0; 0,0,0,0,0];
					[p,~] = size(W);
					alpha = 0;
					lam = 1;
					G = (1:p)';
					groups = (1:p);
					y = soslasso_shrink_logistic(W,G,groups,lam,alpha);
					Wtemp = sign(W).*max(abs(W) - lam,0);
					testCase.assertEqual(norm(y-Wtemp),0,'AbsTol',1e-8);
				end
    end
end
