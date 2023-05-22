function data = imageInterpolator(fileName,xRange,yRange,mode,scale)
%IMAGEINTERPOLATOR allows you to interpolate a graph in an image by
% plotting it in a figure and making you select points on it
%
% Property of THRUST, unauthorized distribution is not allowed
% version: 1.1 - 03/12/2022 - Author: Alessandro Rampazzo
%
%   DATA = IMAGEINTERPOLATOR(FILENAME, XRANGE, YRANGE) outputs a
%   matrix containing a list of the points selected. FILENAME is the name
%   or the path to the image to plot, XRANGE and YRANGE are the ranges of
%   the image (for example the x axis in the image goes from 0 to 1 and the
%   y axis from 1 to 10) given in this format: [low, high]
%
%   DATA = IMAGEINTERPOLATOR(FILENAME, XRANGE, YRANGE, MODE). MODE is a
%   string being "get" (or "g") if you want to select points and or "view" 
%   (or "v") if you want to plot only the image.
%
%   DATA = IMAGEINTERPOLATOR(FILENAME, XRANGE, YRANGE, MODE, SCALE). SCALE
%   is the type of scale of the graph to be plotted, it can be "linear" or
%   "log"
%
%   Example: simple data extraction and plotting
%
%     fileName = "image1.png";
%     xRange = [0,1];
%     yRange = [1,10];
%     % get interpolating data
%     data = imageInterpolator(fileName,xRange,yRange)
%     % plot image
%     imageInterpolator(fileName,xRange,yRange,"v");
%     % plot interpolating data
%     plot(data(:,1),data(:,2))
%
%   See also INTERP1, LINSPACE, POLYFIT, POLYVAL


if ~exist("mode","var")
    mode = "get";
end
if ~exist("scale","var")
    scale = "linear";
end

%% plotting image
f = figure;
hold on

if lower(scale) == "log"
    set(gca,"YScale","log")
    set(gca,"XScale","log")
    ximg = logspace(log10(xRange(1)),log10(xRange(2)));
    yimg = logspace(log10(yRange(1)),log10(yRange(2)));
elseif lower(scale) == "linear"
    ximg = linspace(xRange(1),xRange(2));
    yimg = linspace(yRange(1),yRange(2));
else
    error("scale option not recognized")
end

img = imread(fileName);
img = img(end:-1:1,:,:);
image(ximg,yimg,img)
axis tight

%% abilitating acquisition of data
if lower(mode) == "g" || lower(mode) == "get"
    i = 0;
    data = [];
    p = {};
    set(gcf,"WindowState","fullscreen")

    while true
        [x,y,button] = ginput(1);
        if get(gcf,"WindowState") == "normal"
            break
        end
        if button == 1
            i = i+1;
            data(i,:) = [x,y];
            p{i} = plot(x,y,"ro","markerfacecolor","r");
        elseif button == 3
            i = max(0,i-1);
            if i > 0
                data(end,:) = [];
                p{end}.Visible = "off";
                p = p(1:end-1);
            else
                data = [];
                p{1}.Visible = "off";
                p = [];
            end
        end
    
    end

%     for i = 1:size(data,1)
%         fprintf("%1.8f %1.8f\n",data(i,1),data(i,2))
%     end
    close(f)
end
end