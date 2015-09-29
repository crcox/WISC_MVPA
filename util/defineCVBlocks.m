function cvb = defineCVBlocks(stimcode, varargin)
%DEFINECVBLOCKS  Define k-fold or LOO cross validation assignments
%
% defineCVBlocks(stimcode, 'LOO')
% defineCVBlocks(stimcode, 'folds', k)
% defineCVBlocks(stimcode, 'folds', k, 'schemes', n)
%
% stimcode is either a numeric vector or a cell array of strings that
% identify each stimulus. Sometimes the same stimulus is presented multiple
% times, either over consecutive timepoints, or in separate trials. Exact
% repetitions should be witheld as a set---that is, they should all be
% assigned to the same cross validation block. Leave-one-out cross
% validation will be over stimuli, not trials.
%
% The schemes parameter allow one to specify the number of unique cross
% validation schemes to generate. If n>1, cvb will be a p by n matrix,
% where p is the number of trials.
%
% In an attempt to ensure uniqueness, the algorithm will make n * stoploss
% attempts, where the default stoploss is 10. This just prevents an
% infinite loop.
% 
% Examples:
%   stimcode = repmat(1:8,1,2);
%   defineCVBlocks(stimcode,'folds',4,'schemes',4)
%   % Numeric codes need not be consecutive
%   stimcodeRand = repmat(randperm(100,8),1,2);
%   defineCVBlocks(stimcodeRand,'folds',4,'schemes',4)
%   % Also works with strings
%   stimcodeCellStr = repmat({'a','b','c','d','e','f','g','h'},1,2);
%   defineCVBlocks(stimcodeCellStr,'folds',4,'schemes',4)
%   % Protection against asking for too many unique schemes...
%   defineCVBlocks(stimcodeCellStr,'folds',4,'schemes',10000)
%   % LOO is useful in these more complex cases
%   defineCVBlocks(stimcodeCellStr,'LOO', true)
% See also:
% BAR
% SOMECLASS/SOMEMETHOD

p = inputParser();
addRequired(p, 'stimcode', @isvector);
addParameter(p, 'balance', @isnumeric); % This is not implemented. Should allow balancing over conditions/scanner runs.
addParameter(p, 'LOO', false, @islogical);
addParameter(p, 'folds', 0, @isintegervalue);
addParameter(p, 'schemes', 1, @isintegervalue);
addParameter(p, 'stoploss', 10, @isintegervalue);
parse(p, stimcode, varargin{:});

% Check that stoploss is >0
if p.Results.stoploss < 1
  err.message = sprintf('stoploss must by >0');
  err.identifier = 'WholeBrain_RSA:defineCVBlocks:invalidArgument';
  error(err);
end
stoploss = p.Results.stoploss;

% Check that stimcode is either a numeric vector or a cell array of
% strings. Vector has already been enforced by @isvector.
if ~(isnumeric(p.Results.stimcode) || iscellstr(p.Results.stimcode))
  err.message = sprintf('stimcode must satisfy isnumeric or iscellstr.');
  err.identifier = 'WholeBrain_RSA:defineCVBlocks:invalidArgument';
  error(err);
end

stimcode = p.Results.stimcode;

% Parse stimcode and generate stimid (internal, integer representation).
stimset = unique(stimcode, 'sorted');
stimid = zeros(numel(stimcode),1);
for i = 1:numel(stimset)
  if iscellstr(stimset)
    z = strcmp(stimset{i}, stimcode);
  else %isnumeric
    z = stimset(i) == stimcode;
  end
  stimid(z) = i;
end

% Check that one and only one method is specified.
if p.Results.LOO && (p.Results.folds > 1)
  err.message = sprintf('Cannot specify both LeaveOneOut and kFold.');
  err.identifier = 'WholeBrain_RSA:defineCVBlocks:incompatibleArguments';
  error(err);
elseif ~(p.Results.LOO || (p.Results.folds > 1))
  err.message = sprintf('No method specified.');
  err.identifier = 'WholeBrain_RSA:defineCVBlocks:tooFewArguments';
  error(err);
end

% Check method specifications, and catch invalid cases.
if p.Results.LOO;
  method = 'LeaveOneOut';
elseif p.Results.folds > 1
  method = 'kFold';
  folds  = p.Results.folds;
  % Check that number of schemes is >0
  if p.Results.schemes < 1
    err.message = sprintf('Invalid number of schemes: %d. Must be integer >0.', p.Results.schemes);
    err.identifier = 'WholeBrain_RSA:defineCVBlocks:invalidArgument';
    error(err);
  end
  schemes = p.Results.schemes;
  % Check that number of stimuli can satisfy number of folds.
  if max(stimid) < folds;
    err.message = sprintf('More folds (%d) than unique stimuli (%d).', max(stimid));
    err.identifier = 'WholeBrain_RSA:defineCVBlocks:incompatibleArguments';
    error(err);
  end
else
  err.message = sprintf('Invalid number of folds: %d. Must be integer >1.', p.Results.folds);
  err.identifier = 'WholeBrain_RSA:defineCVBlocks:invalidArgument';
  error(err);
end

n = uint32(max(stimid));
p = uint32(length(stimid));

switch method
  case 'LeaveOneOut'
    % Leave One Out doesn't actually involve anything more than cross
    % validation by stimulus. So just return stimid.
    cvb = stimid;
  case 'kFold'
    % kFold cross validation involves:
    % Determining the size of each fold, accounting for cases where the
    % number of trials is not evenly divisible by the number of folds.
    % Assigning stimuli to folds at random.
    cvb = zeros(schemes, p); % start with this orientation to appease unique(...,'rows');
    m = idivide(p, folds);
    r = mod(p, folds);
    foldsize = ones(1, folds, 'uint32') * m;
    foldsize(1:r) = m + 1;
    scheme = 1;
    iter = 0;
    BAILOUT = false;
    while scheme <= schemes
      if iter > stoploss
        BAILOUT = true;
        break
      end
      
      ix = randperm(n);
      foldix = mat2cell(ix, 1, foldsize);
      for fold = 1:folds
        z = ismember(stimid, foldix{fold});
        cvb(scheme, z) = fold;
      end
      % Check uniqueness
      [~,ia] = unique(cvb(1:scheme,:), 'rows');
      if numel(ia) < scheme
        % roll back, and do not increment scheme.
        cvb(scheme, :) = 0;
      else
        scheme = scheme + 1;
      end
      iter = iter + 1;
    end
    % Transpose to expected orientation.
    cvb = cvb(1:scheme-1,:)';
    if BAILOUT
      warning('Unable to generate %d unique schemes based on input. Returning %d schemes.', schemes, scheme-1);
    end
end
end

% This should be part of a package, probably:
function b = isintegervalue(x)
  b = isa(x,'integer') || (imag(x)==0 && mod(x,1)==0);
end