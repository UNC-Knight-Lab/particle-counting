%% Create histogram of 2nm bins from 0 to 50 nm

X = Combined_measurement_file %insert name of combined measurement file made from Script_Polymer_Analysis
bin_centers = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29,31,33,35,37,39,41,43,45,47,49]; %generate 2 nm bin columns 
edges = [0:2:50]; %make edges between 0 and 50 spaced by intervals of 2
hist = histcounts(X, edges); %function to generate histogram to check particle diatmeter, 

%% Translate row files to columns 

bin_centers_column = bin_centers.' %translates row file to column
hist_column = hist.' %translates row file to column

%% Create normalized histogram text file of measurement file plotting bin centers 

Total_Freq = sum(Combined_Hist_Columns,2);
Freq_max = max(Total_Freq);
Normal_Freq = Total_Freq/Freq_max;
Table = table (bin_centers_column, Normal_Freq);
writetable(Table,'nBA5_XTEN2_sonication_nospaghet_2nmbins.txt','Delimiter','space');
%% Histogram figure 
figure();
edges =[0:2:50];
histogram(X,edges);
title('Put Your Title Here');
set(gca,'box','off');
xlabel('Diameter (nm)');
ylabel('Count');

%% Finding average particle size and standard deviation
Particle_mean = mean(X)
Particle_std = std(X)
