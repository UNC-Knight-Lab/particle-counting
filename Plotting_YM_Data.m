%% Loading in Data File 
cd('D:')    
M = readmatrix('qi-data-2022.12.07-18.23.34.764_processed-2022.12.07-18.26.27.txt');
M(:,1) = [];
M(isnan(M))=0;
%% Exracting XY Coordinates and Young's Modulus
X_Position = M(:,2);
Y_Position = M(:,3);
Youngs_Modulus = M(:,7);
Youngs_Modulus_MPa = Youngs_Modulus/10e6;
%% Finding Area of YM measurement
xmin = min(X_Position);
xmax = max(X_Position);
ymin = min(Y_Position);
ymax = max(Y_Position);
%% Filtering out Mica-Cantilever forces
Youngs_Modulus_Filtered = M(:,7);
Youngs_Modulus_Filtered_MPA = Youngs_Modulus_Filtered/10e6;
Youngs_Modulus_Filtered_MPA(Youngs_Modulus_Filtered_MPA>=100) = 0;

%% Plotting Filtered YM over AFM image 
I = imread('qi-fit-2022.12.07-18.23.34.764_reference-force-height-default.png');
img1 = flip(I);
image(I,'XData',[xmin xmax], 'YData', [ymin ymax]);
hold on
%YM_Data_Filtered = figure();
h = plot3(X_Position,Y_Position,Youngs_Modulus_Filtered_MPA);
title('Filtered YM Data');
set(h, {'LineWidth'}, {0.5});
hold off
%% Histograms of YM Data 
% CHA= CHA_XTEN2_193(:,3);
% edges =(1:2:100);
% g = figure();
% g1 = histogram(CHA,edges);
% title('Particle YM Measurements')

%% Trying to remove zeros from array
% tBA_filter = sort(tBA);
% edges =(1:1:20);
% g = figure();
% g1 = histogram(tBA_filter,edges);
% title('Particle YM Measurements_Filtered')