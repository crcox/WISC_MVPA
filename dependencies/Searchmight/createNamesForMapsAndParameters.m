% create unique names for each map/parameter combination

function [mapNames] = createNamesForMapsAndParameters(maps,mapParameters)

nmaps = length(maps);
mapNames = cell(nmaps,1);

for im = 1:nmaps
  map = maps{im};
  parameters = mapParameters{im};
  
  % convert nice map name to classifier name in the function
  switch map
   case {'voxelwiseGNB','voxelwiseGNBsmoothed','searchlightGNB','searchlightGNB','searchlightGNBsmoothed','searchmightGNB','searchlightLDA','searchlightLDA_ridge','searchlightLDA_shrinkage','searchlightQDA_shrinkage'}
    % map has no parameters
    mapNames{im} = map;
   case {'searchlightSVM_linear','searchlightSVM_quadratic','searchlightSVM_sigmoid','searchlightSVM_rbf'}
    np = length(parameters);
    txt = {};
    
    ip = 1; it = 1;
    while ip <= np
      argval = parameters{ip}; ip = ip + 1;
      switch argval
       case {'lambda'}
        lambda = parameters{ip}; ip = ip + 1;
        switch lambda
         case {'crossValidation'}
          txt{it} = sprintf('lambda-xv'); it = it + 1;
         otherwise
          txt{it} = sprintf('lambda-%s',num2str(lambda)); it = it + 1;
        end
       case {'gamma'}
        gamma = parameters{ip}; ip = ip + 1;
        switch lambda
         case {'crossValidation'}
          txt{it} = sprintf('gamma-xv'); it = it + 1;
         otherwise
          txt{it} = sprintf('gamma-%s',num2str(gamma)); it = it + 1;
        end
       otherwise
        fprintf('error: parameter %s does not go with map %s\n',argval,map);return;
      end
    end
    
    if isempty(txt)
      mapNames{im} = sprintf('%s-defaults',map);
    else
      mapNames{im} = sprintf('%s-%s',map,txt{1});
      for it = 2:length(txt)
        mapNames{im} = sprintf('%s_%s',mapNames{im},txt{it});
      end
    end
   
   case {'searchlightLR_L1','searchlightLR_L2'}
    mapNames{im} = map; % not using these yet
   
   case {'knn'}
    np = length(parameters);
    txt = {};
    
    ip = 1; it = 1;
    while ip <= np
      argval = parameters{ip}; ip = ip + 1;
      switch argval
       case {'k'}
        K = parameters{ip}; ip = ip + 1;
        txt{it} = sprintf('K-%s',num2str(K)); it = it + 1;
       otherwise
        fprintf('error: parameter %s does not go with map %s\n',argval,map);return;
      end
    end
    
    if isempty(txt)
      mapNames{im} = sprintf('%s-defaults',map);
    else
      mapNames{im} = sprintf('%s-%s',map,txt{1});
      for it = 2:length(txt)
        mapNames{im} = sprintf('%s_%s',mapNames{im},txt{it});
      end
    end
    
  end; % switch on classifier
end; % for on map cell array
