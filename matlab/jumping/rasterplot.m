function handle = rasterplot(time,x,width,height)
%RASTERPLOT plots a spike raster plot for a given input
% 
if nargin < 3
    width = 1;
end
if nargin < 4
    height = 1;
end
time =reshape(time,1,[]);
x =reshape(x,1,[]);
N = length(x);

xx = [x;x+height;nan(1,N)];
time = [time;time;time];

xx = xx(:);
time = time(:);

handle = line(time,xx,'LineWidth',width);
end