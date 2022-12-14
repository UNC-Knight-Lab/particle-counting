%% 
bin_centers = [1,3,5,7,9,11,13,15,17,19,21,23,25,27,29]; %generate 1 column vector for bin centers
edges = (0:2:30);
% tBAXTEN2_10 = table2array(tBAXTEN210measurementfile);
% tBAXTEN2_12 = table2array(tBAXTEN212measurementfile);
% tBAXTEN2_18 = table2array(tBAXTEN218measurementfile);
% tBAXTEN2_19 = table2array(tBAXTEN219measurementfile);
% tBAXTEN2_9 = table2array(tBAXTEN29measurementfile);

%% 

bin_centers_column = bin_centers.'
% CHAXTEN2_10_Frequency = histcounts(CHA_XTEN2_10_measurementfile,edges);
% CHAXTEN2_10_column = CHAXTEN2_10_Frequency.';
% CHA_XTEN2_12_Frequency = histcounts(CHA_XTEN2_12_measurementfile,edges);
% CHA_XTEN2_12_column = CHA_XTEN2_12_Frequency.';
% CHA_XTEN2_7_Frequency = histcounts(CHA_XTEN2_7_measurementfile,edges);
% CHA_XTEN2_7_column = CHA_XTEN2_7_Frequency.';
EA_XTEN2_17_Frequency = histcounts(EA_XTEN2_17_measurementfile,edges);
EA_XTEN2_17_column = EA_XTEN2_17_Frequency.';
%%                                                                                                                                                                     
Combined_Hist_Columns = [EA_XTEN2_17_column];
Total_Freq = sum(Combined_Hist_Columns,2);
Table = table (bin_centers_column, Total_Freq);
writetable(Table,'EA-XTEN2.txt','Delimiter','space')
%% 


%figure()
%edges =[20:2:50]
%histogram(Freq_1321_column,edges)


%% Finding average particle size and standard deviation
