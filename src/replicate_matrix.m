function [A_tilde, groups, group_arr] = replicate_matrix(A,group_MAT)

% Function to form 'A_tilde' that is a replicated version of A
% for the replication strategy, please refer to
% Jacob et. al. Group Lasso with Overlaps and Graph lasso, ICML 2009

% INPUTS:
% A         = Original Data
% group_Mat = cell array having rows as groups


% OUTPUTS:
% A_tilde   = replicated data
% groups = indices of data and group to which it belongs
% group_arr = matrix whose rows correspond to group membership

  if ~iscell(A)
    l = length(group_MAT);

    % REPLICATION STEP
    idx = 1:l;

    groups = [];
    A_tilde = [];
    group_start = zeros(l,1); group_len = zeros(l,1);
    group_start(1)=1;
    % groups replication, and A_tilde if replication used (by replicating columns of A)
    for i = idx

      subg = repmat(i,1,length(group_MAT{i}));
      %replicate columns of A if forcing not used
      subc = A(:,group_MAT{i});
      A_tilde = [A_tilde subc];

      groups = [groups subg];
      if (i>1)
        group_start(i) = group_start(i-1) + group_len(i-1);
      end
      group_len(i)=length(group_MAT{i});

    end

    % construct group_arr structure, for which row i consists of the
    % indices in group i, padded out with "n+1", which points to a
    % dummy index in x
    group_arr=(length(groups)+1)*ones(l,max(group_len));
    for i=1:l
      group_arr(i,1:group_len(i))=group_start(i):group_start(i)+group_len(i)-1;
    end

  else

    X = A;
    J = length(X);
    A_t = cell(J,1);
    for j = 1:J

      A = X{j};
      l = length(group_MAT);

      % REPLICATION STEP
      idx = 1:l;

      %		 groups = [];
      %		 A_tilde = [];

      %preallocate A_tilde
      if j==1
        group_start = zeros(l,1); group_len = zeros(l,1);
        group_start(1)=1;
        for i = idx
          group_len(i)=length(group_MAT{i});
        end

      end
      A_tilde = zeros(k,sum(group_len));
      if j == 1
        groups  = zeros(sum(group_len),1);
      end


      % groups replication, and A_tilde if replication used (by replicating columns of A)
      for i = idx


        subg = repmat(i,1,group_len(i));
        %replicate columns of A if forcing not used
        subc = A(:,group_MAT{i});
        %			 A_tilde = [A_tilde subc];
        %			 groups = [groups subg];
        if j==1
          if (i>1)
            group_start(i) = group_start(i-1) + group_len(i-1);
          end
        end
        A_tilde(:,group_start(i):group_start(i)-1+group_len(i)) = subc;
        groups(group_start(i):group_start(i)-1+group_len(i)) = subg;

      end

      % construct group_arr structure, for which row i consists of the
      % indices in group i, padded out with "n+1", which points to a
      % dummy index in x
      if j==1
        group_arr=(length(groups)+1)*ones(l,max(group_len));
        for i=1:l
          group_arr(i,1:group_len(i))=group_start(i):group_start(i)+group_len(i)-1;
        end
      end

      A_t{j} = A_tilde;

    end

    A_tilde = A_t;
  end

end
