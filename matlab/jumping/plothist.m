% lap detector
% SGL 2020-11-28
clc; close all;
if ~exist('pos','var')
    exp_directory = '~/onedrive/JHU/913_Jumping_Recording/2020-11-11_Rat913-02';
    mat_filename = fullfile(exp_directory,'data.mat');
    load(mat_filename,'pos', 'spike')
end
% plot(pos.x, pos.y,'.')
% axis equal

figure(2)
% recording
r1i = 105.221068;
r1f = 703.846944;
r2i = 945.026778;
r2f = 1341.248894;

xmax = 166;
rectangle('Position',[r1i,0,r1f-r1i,xmax],'FaceColor','y')
hold on
rectangle('Position',[r2i,0,r2f-r2i,xmax],'FaceColor','y')

% maze
m1i = 408.650168;
m1f = 677.921799;
m2i = 988.727761;
m2f = 1303.765711;
   
rectangle('Position',[m1i,0,m1f-m1i,xmax],'FaceColor',[0 .8 .8])
hold on
rectangle('Position',[m2i,0,m2f-m2i,xmax],'FaceColor',[0 .8 .8])

% lap detection range
t1i = 110;
t1f = 659;
t2i = 945;
t2f = 1282;

t = pos.t(((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f)));
x = pos.x(((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f)));
y = pos.y(((pos.t > t1i) & (pos.t < t1f)) | ((pos.t > t2i) & (pos.t < t2f)));

plot(t, x, '.')

x_thresh = 82;

%plot(t, xmax * (x > x_thresh), '.-r')

ylim([0 xmax])

% times that the rats jump (left to right and right to left)
j1 = t(diff(x > x_thresh) == -1);
%j2 = t(diff(x > x_thresh) == 1);

N = length(j1)-1;
lr = ones(N,1); % lap to right
ll = ones(N,1); % lap to left
x1 = ones(N,1);
x2 = ones(N,1);
    
for i=1:N
    % lap
    tl = t(t>=j1(i) & t<=j1(i+1));
    [x1(i), idx] = min(x(t>j1(i) & t<j1(i+1)));
    lr(i) = tl(idx);
    [x2(i), idx] = max(x(t>j1(i) & t<j1(i+1)));
    ll(i) = tl(idx);
end

% exception for beginning and end
tl = t(t<j1(1));
[x20, idx] = max(x(t<j1(1)));
l20 = tl(idx);
x2 = [x20;x2];
ll = [l20;ll];
tl = t(t>j1(N+1));
[x1(N+1), idx] = min(x(t>j1(N+1)));
lr(N+1) = tl(idx);

plot(lr, x1, 'pk', 'MarkerSize',15)
plot(ll, x2, 'sk', 'MarkerSize',15)

plot(j1, x_thresh * ones(size(j1)), 'hk', 'MarkerSize',15)

figure(3)
for i = 1:N+1
   subplot(N+1,2,2*i-1);
   plot(x, t)
   xlim([0 xmax])
   ylim([ll(i) lr(i)])
   if i==N+1
       break;
   end
   subplot(N+1,2,2*i);
   plot(x, t)
   xlim([0 xmax])
   ylim([lr(i) ll(i+1)])
end