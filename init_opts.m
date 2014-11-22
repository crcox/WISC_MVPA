%% FUNCTION init_opts
% initialization options for multi-task learning library
%
% If one of the ncessary opts are empty then it will be set to default
% values as specified in this file.
%
% Table of Options.  * * indicates default value.
% FIELD                DESCRIPTION
%% Optimization options
%
%  .max_iter               Maximum iteration step number
%                           *1000*
%  .tol                    Tolerance
%                           *10e-3*
%  .tFlag                  Termination condition
%                           0 => change of absolute function value:
%                             abs( funcVal(i)- funcVal(i-1) ) <= .tol
%                         * 1 => change of relative function value:
%                             abs( funcVal(i)- funcVal(i-1) )
%                              <= .tol * funcVal(i-1)
%                           2 => absolute function value:
%                             funcVal(end)<= .tol
%                           3 => Run the code for .maxIter iterations
%% Starting Point
%
% .W0               Starting point of W.
%                   Initialized according to .init.
%
% .C0               Starting point for the intercept C (for Logistic Loss)
%                   Initialized according to .init.
%
% .init             .init specifies how to initialize W, C.
%                         1 => .W0 and .C0 are defined
%                       * 2 => .W0= zeros(.), .C0=0 *
%
%
%% LICENSE
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%   Copyright (C) 2011 - 2012 Jiayu Zhou and Jieping Ye
%
%   You are suggested to first read the Manual.
%   For any problem, please contact with Jiayu Zhou via jiayu.zhou@asu.edu
%
%   Last modified on June 21, 2012.
%

function opts = init_opts (opts)

%% Default values
DEFAULT.max_iter   = 1000;
DEFAULT.tol        = 1e-8;
MINIMUM_TOLERANCE  = eps * 100;
DEFAULT.tFlag      = 1;
DEFAULT.init       = 2;
DEFAULT.Debias     = true;
DEFAULT.GroupSize  = 1;
DEFAULT.GroupShift = DEFAULT.GroupSize;
DEFAULT.SharedSpaceVoxelSize = 1; % must be isotropic.
DEFAULT.alpha      = [];
DEFAULT.lambda     = [];

%% Starting Point
curField = 'init';
if isfield(opts,curField)
    if (opts.(curField)~=0) && (opts.(curField)~=1) && (opts.(curField)~=2)
        opts.(curField)=DEFAULT.(curField); % if .init is not 1, or 2, then use the default.
    end
    
    if (~isfield(opts,'W0')) && (opts.init==1)
        opts.(curField)=DEFAULT.(curField); % if .W0 is not defined and .init=1, set .init=default
    end
else
    opts.(curField) = DEFAULT.(curField); % if .init is not specified, use default
end

%% Tolerance
curField = 'tol';
if isfield(opts, curField)
    % detect if the tolerance is smaller than minimum
    % tolerance allowed.
    if (opts.(curField) < MINIMUM_TOLERANCE)
        opts.(curField) = MINIMUM_TOLERANCE;
    end
else
    opts.(curField) = DEFAULT.(curField);
end

%% Maximum iteration steps
curField = 'max_iter';
if isfield(opts, curField)
    if (opts.(curField)<1)
        opts.(curField) = DEFAULT.(curField);
    end
else
    opts.(curField) = DEFAULT.(curField);
end

%% Termination condition
curField = 'tFlag';
if isfield(opts,curField)
    if (opts.tFlag~=1) && (opts.tFlag~=2)
        opts.tFlag=DEFAULT.(curField);
        IllegalOptWarning(opts.(curField), curField, DEFAULT.(curField));
    else
        opts.tFlag=floor(opts.tFlag);
    end
else
    opts.tFlag = DEFAULT.(curField);
end

%% Grouping Options
curField = 'GroupSize';
if isfield(opts,curField)
    opts.(curField) = uint32(opts.(curField));
    if any(opts.(curField)==0)
        opts.(curField)=DEFAULT.(curField);
        IllegalOptWarning(opts.(curField), curField, DEFAULT.(curField));
    end
    if length(opts.(curField)) == 1;
        opts.(curField)=repmat(opts.(curField),1,3);
    elseif length(opts.(curField)) > 3;
        error('Specified >3 dimensions for opts.GroupSize.')
    end
else
    opts.(curField) = DEFAULT.(curField);
end

curField = 'GroupShift';
if isfield(opts,curField)
    opts.(curField) = opts.(curField);
    if any(opts.(curField)==0)
        opts.(curField)=DEFAULT.(curField);
        IllegalOptWarning(opts.(curField), curField, DEFAULT.(curField));
    end
    if length(opts.(curField)) == 1;
        opts.(curField)=repmat(opts.(curField),1,3);
    elseif length(opts.(curField)) > 3;
        error('Specified >3 dimensions for opts.GroupShift.')
    end
else
    opts.(curField) = DEFAULT.(curField);
end

curField = 'SharedSpaceVoxelSize';
if isfield(opts,curField)
    opts.(curField) = opts.(curField);
    if any(opts.(curField)==0)
        opts.(curField)=DEFAULT.(curField);
        IllegalOptWarning(opts.(curField), curField, DEFAULT.(curField));
    end
% %     THIS CONSTRAINT IS NO LONGER NEEDED
%     if length(opts.(curField)) > 1;
%         error('Specified >1 value for opts.SharedSpaceVoxelSize.')
%     end
else
    opts.(curField) = DEFAULT.(curField);
end

%% Regularizer parameters
curField = 'alpha';
if isfield(opts,curField)
    if any([opts.(curField)<0,opts.(curField)>1])
        opts.(curField)=DEFAULT.(curField);
        error('opts.alpha must be a vector of scalars in the range 0 to 1; endpoints correspond to lasso and group lasso, respectively.');
    end
else
    opts.(curField) = DEFAULT.(curField);
    warning('remember to set opts.alpha values between 0 and 1.')
end

curField = 'lambda';
if isfield(opts,curField)
    if opts.(curField)<0
        opts.(curField)=DEFAULT.(curField);
        error('opts.lambda must be a vector of scalars >0.')
    end
else
    opts.(curField) = DEFAULT.(curField);
    warning('remember to set opts.lambda values >0.')
end

%% Misc
curField = 'Debias';
if isfield(opts,curField)
    if (opts.(curField)~=0) && (opts.(curField)~=1)
        opts.(curField)=DEFAULT.(curField);
        IllegalOptWarning(opts.(curField), curField, DEFAULT.(curField));
    else
        opts.(curField)=logical(opts.(curField));
    end
else
    opts.(curField) = DEFAULT.(curField);
end

end


function IllegalOptWarning(val,field,dval)
    warning('Illegal value %d for opts.%s... assigning default (%d).', ...
        val, field, dval);
end
