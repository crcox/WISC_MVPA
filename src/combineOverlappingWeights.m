function W = combineOverlappingWeights(Wc, G)
  % Create indexes into the replicated weight matrix.
  S = subjectfilters(G);
  % Combine (and, in the process, sort into subject order).
  N = cellmax(G);
  W = cell(1,size(G,2));
  fprintf('%8s%8s%8s\n','subject','nnz','unique');
  for j = 1:size(G,2)
    W{j} = zeros(N(j),1);
    for i = 1:size(G,1); 
      g = G{i,j};
      s = S{i,j};
%     s = group_arr(k,:)+1; % to account for bias weight
      W{j}(g) = W{j}(g) + Wc(s,j);
    end
    s = cell2mat(S(:,j));
    w = Wc(s,j);
    t = cell2mat(G(:,j));
    fprintf('% 8d% 8d% 8d\n', j, nnz(w), numel(unique(t(w~=0))));
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
  N = uint32(max(cellfun(@max, C)));
end
