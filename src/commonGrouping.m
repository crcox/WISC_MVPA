function [Gc,ix] = commonGrouping(G)
    n = uint32(max(cellfun('length', G),[],2));
    c = uint32([0;cumsum(n(:))]);
    N = uint32(max(cellfun(@max, G)));
    for i = 1:size(G,1)
        for j = 1:size(G,2)
            x = uint32(G{i,j});
            g = uint32([x; ones(n(i)-length(x),1,'uint32')*(N(j)+1)]);
            G{i,j} = g;
        end
    end
    ix = cell2mat(G);
    
    Gc = cell(size(G,1),1);
    for i = 1:length(Gc)
        a = c(i) + 1;
        b = c(i + 1);
        Gc{i} = uint32(a:b);
    end
end
    
%     N = max(cellfun(@max, G));
%     for i = 1:size(G,1)
%         for j = 1:size(G,2)
%             x = G{i,j};
%             G{i,j} = [x; ones(n(i)-length(x),1)*(N(j)+1)];
%         end
%     end
%     
%     ix = cell2mat(G);
%     z = false(size(ix));
%     for i = 1:size(ix, 2)
%         x = ix(:,i);
%         [~,ai] = unique(x,'stable');
%         z(ai,i) = true;
%         z(x==(N(i)+1),i) = true;
%     end
%     sum(z)
%     ix = reshape(ix(z), [], size(ix, 2)); 
    