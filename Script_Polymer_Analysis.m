%% Pick images for input

maskdir = 'Mask of EA5-XTEN2_17.tif'; %this is the black and white mask
cntmaskdir = 'Count Masks of EA5-XTEN2_17.tif'; %this is the countmask
rawim = 'EA5-XTEN2_17.tif'; %this is the raw image
pixelsize = 6.48;


tmin = 1;
tmax =[];
tmax = size(imfinfo(maskdir), 1);



%% load images in
for c = 1:3;
    if c == 1
        gadir = maskdir;
        ctdir = cntmaskdir;
    end  
         
    for i = 1:tmax;
        preimin{i,1} = imread(gadir);
        preiminct{i,1} = imread(ctdir);
        preiminraw{i,1} = imread(rawim);
    end

    for i = 1:tmax
        if c ==1;
            maskin{i,1} = preimin{i,1};
            maskinct{i,1} = preiminct{i,1};
        end    
    end
end

figure();
imagesc(maskin{1,1});
figure();
imagesc(maskinct{1,1});

%% Gather some data about the images
tnum = tmax - tmin + 1;
max_field = zeros(tnum,1);
peak_thresh = cell(tmax,1);
currmax = [];
low_thresh_clean = cell(tnum,1);
low_thresh_count = zeros(tnum,1);
peak_thresh_count = zeros(tnum,1);
labeledmask_peak = cell(tnum,1);
labeledmask = cell(tnum,1);


%% Count number of objects in mask

for i = 1:tmax;
    labeledmask{i,1} = bwlabel(maskin{i,1});
    cellcount(i,1) = max(max(labeledmask{i,1}));
end

figure()
imagesc(labeledmask{1,1});

%% Find the indices of every pixel in the object, which is an indication of area

total = 1;
for i = 1:tnum;
    for j = 1:cellcount(i,1);
        currmask = [];
        currmask = labeledmask{i,1}==j;
        [xs,ys] = find(currmask); 
        ind{i,j} = [ys,xs];
    end
end


%% Fit an elipse to each object and find the major and minor diameter

measurements = regionprops(labeledmask{1,1},{'Centroid','Orientation','MajorAxisLength','MinorAxisLength'});

%% Extract data to matrices for the centroid, major axis, and minor axis of the ellipse
majax =[];
minax =[];
majaxnm =[];
minaxnm =[];

data = struct2cell(measurements);

for i = 1:size(data,2)
    majax(i,1) = data{2,i};
    minax(i,1) = data{3,i};
    centroidxy(i,1) = data{1,i}(1,1);
    centroidxy(i,2) = data{1,i}(1,2);
end

majaxnm = majax./pixelsize;
minaxnm = minax./pixelsize;

EA_XTEN2_17_measurementfile = sort(majaxnm)
%figure()
%edges =[0:2:30]
%histogram(majaxnm_1,edges)
%hold on
%histogram(minaxnm,edges)

%writematrix(majaxnm,'major_axis_diameter_1014.txt','Delimiter','tab')
%type 'major_axis_diameter.txt'

%writematrix(minaxnm,'minor_axis_diameter_1014.txt','Delimiter','tab')
%type 'minor_axis_diameter.txt'

%% Plot each ellpise on the mask and on the original image


figure()
imagesc(labeledmask{1,1})
hold on
t = linspace(0,2*pi,50);
for l = 1:length(measurements)
    a = measurements(l).MajorAxisLength/2;
    b = measurements(l).MinorAxisLength/2;
    Xc = measurements(l).Centroid(1);
    Yc = measurements(l).Centroid(2);
    phi = deg2rad(-measurements(l).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'r','Linewidth',1)
end
hold off


fig =figure;
imagesc(preiminraw{1,1})
hold on
t = linspace(0,2*pi,50);
for l = 1:length(measurements)
    a = measurements(l).MajorAxisLength/2;
    b = measurements(l).MinorAxisLength/2;
    Xc = measurements(l).Centroid(1);
    Yc = measurements(l).Centroid(2);
    phi = deg2rad(-measurements(l).Orientation);
    x = Xc + a*cos(t)*cos(phi) - b*sin(t)*sin(phi);
    y = Yc + a*cos(t)*sin(phi) + b*sin(t)*cos(phi);
    plot(x,y,'y','Linewidth',1)
end

%% %% Labeling Particle Counting with Original Image  

%Converting the centroid coordinates into an X and Y vector 
Xcd = centroidxy(:,1)
Ycd = centroidxy(:,2)

%Loading the original image and labeling particles 
fig =figure;
imagesc(preiminraw{1,1})
hold on
for ii = 1:length(Ycd)
    text(Xcd(ii),Ycd(ii),num2str(ii),'Color','r')
end
plot(Xcd,Ycd,'o') 
%% Making Histograms of Particle Measurements 

h = histogram(majaxnm);
%xlabel('Diameter (nm)')
%ylabel('Frequency')

counts = h.Values
bincenters = h.BinEdges + h.BinWidth/2

counts_column = counts.'
bincenters_columns = bincenters.'

Histogram_Table = table(bincenters_columns, counts_column)
%writetable(h,'Hiram.txt','Delimiter','space')

% Obtain the handle of the figure that contains the histogram
handleOfHistogramFigure = ancestor(h, 'figure');
% Make the figure window visible in case it was invisible before
handleOfHistogramFigure.Visible  = 'on'
% Bring the figure window to the front
figure(handleOfHistogramFigure);


%% 
edges =(1:2:30);
g = figure();
g1 = histogram(EA_XTEN2_17_measurementfile,edges);


