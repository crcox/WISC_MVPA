ii = 0;
figure(2)
lambdatick = cell(6,1);
lambdalab = cell(6,1);
[lambdatick{[1,6]}] = deal(linspace(0,100,11));
[lambdatick{2:5}] = deal(linspace(0,60,11));
[lambdalab{[1,6]}] = deal(arrayfun(@(x) sprintf('%.0f', x), linspace(0,100,11), 'Unif', 0));
[lambdalab{2:5}] = deal(arrayfun(@(x) sprintf('%.0f', x), linspace(0,60,11), 'Unif', 0));
alpha = linspace(0,1,6);
for i = [1,3,4,5,6,7];
  Bz2 = OLD{i};
  Yz = cell(1,numel(Bz2));
  for j = 1:numel(Bz2)
    Yz{j} = zeros(nitem,nsubj);
    bz = Bz2{j};
    for i = 1:nsubj
      Yz{j}(:,i) = X{i} * bz(:,i);
    end
  end
  ErrorC = cellfun(errc, Yz, 'Unif', 0);
  ii = ii + 1;
  subplot(2,3,ii);
  plot(cell2mat(ErrorC'))
  title(sprintf('%.1f',alpha(ii)));
  set(gca, 'XTick', (lambdatick{ii}/max(lambdatick{ii}))*100);
  set(gca, 'XTickLabel', lambdalab{ii});
  xlabel('lambda')
  ylabel('error');
end
plot(cell2mat(ErrorC'))
ToPlot;
