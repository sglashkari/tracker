
close all
ddd = randi(20,5);
ddd = sort(ddd(:)/25);
figure(1); histogram(ddd,'BinWidth',0.0333);
[counts5,edges] = histcounts(ddd,'BinWidth',0.0333);
hold on; plot(edges(2:end),counts5,'-or')
counts5'
edges'
n = 10;
edges50 = movmean(edges,n+1)'
counts50 = movsum(counts5,n)'

fl = floor(n/2);
figure(2); bar(edges50(6:end-5),counts50(5:end-5))
%bar(edges50(4:end-2),counts50(1:end-2))
whos counts5 counts50 edges edges50
ddd
