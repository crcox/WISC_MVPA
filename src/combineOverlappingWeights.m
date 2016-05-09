function W = combineOverlappingWeights(Wc, G, varargin)
  p = inputParser;
  addRequired(p, 'Wc');
  addRequired(p, 'G');
  addParameter(p , 'verbose' , true);
  parse(p, Wc,G,varargin{:});

  Wc = p.Results.Wc;
  G = p.Results.G;
  verbose = p.Results.verbose;
  % Create indexes into the replicated weight matrix.
  S = subjectfilters(G);
  % Combine (and, in the process, sort into subject order).
  N = cellmax(G);
  W = cell(1,size(G,2));
  if verbose
    fprintf('%8s%8s%8s\n','subject','nnz','unique');
  end
  for j = 1:size(G,2)
    W{j} = zeros(N(j),1);
    for i = 1:size(G,1);
      g = G{i,j};
      s = S{i,j};
%     s = group_arr(k,:)+1; % to account for bias weight
      if ~isempty(g) && ~isempty(s)
        W{j}(g) = W{j}(g) + Wc(s,j);
      end
    end
    s = cell2mat(S(:,j));
    w = Wc(s,j);
    t = cell2mat(G(:,j));
    if verbose
      fprintf('% 8d% 8d% 8d\n', j, nnz(w), numel(unique(t(w~=0))));
    end
	end
end

function S = subjectfilters(G)
  n = cellfun('length', G);
  m = max(n,[],2);
  mc = [0; cumsum(m)];
  S = cell(size(G));
  for j = 1:size(G,2);
    for i = 1:size(G,1);
      a = mc(i) + 1;
      b = mc(i) + n(i,j);
      S{i,j} = (a:b)';
    end
  end
end

function N = cellmax(C)
  C(cellfun('isempty',C)) = {uint32(0)}; %handle empty cells
  N = uint32(max(cellfun(@max, C),[],1));
end
