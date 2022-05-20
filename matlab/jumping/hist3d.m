clc; close all

x =  posi.p(:,1);
y =  posi.p(:,2);
z =  posi.p(:,3);

c = 18;

[~,Xedges] = histcounts(x,'BinWidth',3);
[~,Yedges] = histcounts(y,'BinWidth',3);
[~,Zedges] = histcounts(z,'BinWidth',3);

N = zeros(length(Xedges)-1,length(Yedges)-1,length(Zedges)-1);
for i=1:length(Yedges)-1
    for j=1:length(Zedges)-1
        idx = posi.p(:,2) >= Yedges(i) & posi.p(:,2) < Yedges(i+1) & posi.p(:,3) >= Zedges(j) & posi.p(:,3) < Zedges(j+1);
        N(:,i,j) = histcounts(posi.p(idx,1),Xedges)';
    end
end

X = movmean(Xedges,2); X = X(2:end);
Y = movmean(Yedges,2); Y = Y(2:end);
Z = movmean(Zedges,2); Z = Z(2:end);

[X, Y, Z] = meshgrid(Y, X, Z);

figure
scatter3(X(:),Y(:),Z(:),N(:)/10+eps,'MarkerFaceColor', [0.7 0 0.7])
view([70 14])
axis equal

M = zeros(length(Xedges)-1,length(Yedges)-1,length(Zedges)-1);
for i=1:length(Yedges)-1
    for j=1:length(Zedges)-1
        idx = cluster(c).p(:,2) >= Yedges(i) & cluster(c).p(:,2) < Yedges(i+1) & cluster(c).p(:,3) >= Zedges(j) & cluster(c).p(:,3) < Zedges(j+1);
        M(:,i,j) = histcounts(cluster(c).p(idx,1),Xedges)';
    end
end

figure
scatter3(X(:),Y(:),Z(:),M(:)*10+eps,'MarkerFaceColor',[0.4 0.4 0])
view([70 14])
axis equal

%%
P = (M>max([5 0.2*max(M(:))])).* 5e4.*M./(N+eps);
figure
scatter3(X(:),Y(:),Z(:),P(:)+eps,'MarkerFaceColor',colors(c))
view([70 14])
axis equal