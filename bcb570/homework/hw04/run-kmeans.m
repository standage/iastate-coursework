x = dlmread('wdbc-values.data', ',');
c = kmeans(x, 2);
dlmwrite('wdbc-clusters.data', c);
exit
