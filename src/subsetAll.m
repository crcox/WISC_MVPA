function X = subsetAll(X, rows, cols)
  p = inputParser;
  addRequired(p, 'X');
  addOptional(p, 'rows', [], @islogicallikeOrEmpty);
  addOptional(p, 'cols', [], @islogicallikeOrEmpty);
  parse(p, X, rows, cols);

  X = p.Results.X;
  N = length(X);

  ROWS = processfilter(p.Results.rows,1);
  COLS = processfilter(p.Results.cols,2);
  
  if iscell(X) && ~iscell(ROWS) && ~iscell(COLS)
    X = cellfun(@(x) x(ROWS,COLS), X, 'UniformOutput', false);
  elseif iscell(X) && (iscell(ROWS) || iscell(COLS))
    for i = 1:N
        X{i} = X{i}(ROWS{i}, COLS{i});
    end
  elseif ~iscell(X)
    X = X(ROWS,COLS);
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Process Filter
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function c = processfilter(x, margin)
    if iscell(x)
      assert(N==length(x));
      for ii = 1:N
        sameLength = size(X{ii},margin)==numel(x{ii});
        assert(sameLength);
      end
      c = x;
    else
      if isempty(x)
        c = cell(N,1);
        for ii = 1:N
          c{ii} = true(size(X{ii},margin),1);
        end
      else
        c = repmat({x},N,1);
      end
    end
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function b = islogicallikeOrEmpty(x)
  if iscell(x)
    N = numel(x);
    b = false(N,1);
    for i = 1:N
      b(i) = all(any(bsxfun(@eq, x{i},[1,0]),2));
    end
  elseif isempty(x)
    b = true;
  else
    b = all(any(bsxfun(@eq, x,[1,0]),2));
  end
  b = all(b);
end