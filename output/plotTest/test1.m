clc
close all
clear all

% bands 27282
file_10years = 'bands_51.65625_-116.59375_10years';
file_10yearsWithState = 'bands_51.65625_-116.59375_10yearsWithState';
file_20years = 'bands_51.65625_-116.59375_20years';

% % cell 27282
% file_10years = 'cell_51.65625_-116.59375_10years';
% file_10yearsWithState = 'cell_51.65625_-116.59375_10yearsWithState';
% file_20years = 'cell_51.65625_-116.59375_20years';

% % bands 27283
% file_10years = 'bands_51.65625_-116.53125_10years';
% file_10yearsWithState = 'bands_51.65625_-116.53125_10yearsWithState';
% file_20years = 'bands_51.65625_-116.53125_20years';


% % cell 27283
% file_10years = 'cell_51.65625_-116.53125_10years';
% file_10yearsWithState = 'cell_51.65625_-116.53125_10yearsWithState';
% file_20years = 'cell_51.65625_-116.53125_20years';



% The number of columns in the output file
if file_10years(1) == 'b' % band file is read
%     header = {'YEAR'; 'MONTH'; 'DAY'; 'OUT_SWE_BAND_0'; 'OUT_SWE_BAND_1'; 'OUT_SWE_BAND_2';	'OUT_SWE_BAND_3'; 'OUT_SWE_BAND_4'; 'OUT_SWE_BAND_5'; 'OUT_SWE_BAND_6';	'OUT_SWE_BAND_7'; 'OUT_SWE_BAND_8'; 'OUT_SWE_BAND_9'; 'OUT_SWE_BAND_10'; 'OUT_SWE_BAND_11'; 'OUT_SWE_BAND_12'; 'OUT_GLAC_ACCUM_BAND_0'; 'OUT_GLAC_ACCUM_BAND_1'; 'OUT_GLAC_ACCUM_BAND_2'; 'OUT_GLAC_ACCUM_BAND_3'; 'OUT_GLAC_ACCUM_BAND_4'; 'OUT_GLAC_ACCUM_BAND_5'; 'OUT_GLAC_ACCUM_BAND_6'; 'OUT_GLAC_ACCUM_BAND_7'; 'OUT_GLAC_ACCUM_BAND_8'; 'OUT_GLAC_ACCUM_BAND_9'; 'OUT_GLAC_ACCUM_BAND_10'; 'OUT_GLAC_ACCUM_BAND_11'; 'OUT_GLAC_ACCUM_BAND_12'; 'OUT_GLAC_MELT_BAND_0'; 'OUT_GLAC_MELT_BAND_1'; 'OUT_GLAC_MELT_BAND_2'; 'OUT_GLAC_MELT_BAND_3'; 'OUT_GLAC_MELT_BAND_4'; 'OUT_GLAC_MELT_BAND_5'; 'OUT_GLAC_MELT_BAND_6'; 'OUT_GLAC_MELT_BAND_7';	'OUT_GLAC_MELT_BAND_8';	'OUT_GLAC_MELT_BAND_9'; 'OUT_GLAC_MELT_BAND_10'; 'OUT_GLAC_MELT_BAND_11'; 'OUT_GLAC_MELT_BAND_12'; 'OUT_GLAC_OUTFLOW_BAND_0'; 'OUT_GLAC_OUTFLOW_BAND_1'; 'OUT_GLAC_OUTFLOW_BAND_2'; 'OUT_GLAC_OUTFLOW_BAND_3'; 'OUT_GLAC_OUTFLOW_BAND_4'; 'OUT_GLAC_OUTFLOW_BAND_5'; 'OUT_GLAC_OUTFLOW_BAND_6'; 'OUT_GLAC_OUTFLOW_BAND_7'; 'OUT_GLAC_OUTFLOW_BAND_8'; 'OUT_GLAC_OUTFLOW_BAND_9'; 'OUT_GLAC_OUTFLOW_BAND_10'; 'OUT_GLAC_OUTFLOW_BAND_11'; 'OUT_GLAC_OUTFLOW_BAND_12'; 'OUT_GLAC_MBAL_BAND_0'; 'OUT_GLAC_MBAL_BAND_1'; 'OUT_GLAC_MBAL_BAND_2'; 'OUT_GLAC_MBAL_BAND_3'; 'OUT_GLAC_MBAL_BAND_4'; 'OUT_GLAC_MBAL_BAND_5'; 'OUT_GLAC_MBAL_BAND_6'; 'OUT_GLAC_MBAL_BAND_7'; 'OUT_GLAC_MBAL_BAND_8'; 'OUT_GLAC_MBAL_BAND_9'; 'OUT_GLAC_MBAL_BAND_10'; 'OUT_GLAC_MBAL_BAND_11'; 'OUT_GLAC_MBAL_BAND_12'; 'OUT_GLAC_IMBAL_BAND_0'; 'OUT_GLAC_IMBAL_BAND_1';	'OUT_GLAC_IMBAL_BAND_2'; 'OUT_GLAC_IMBAL_BAND_3'; 'OUT_GLAC_IMBAL_BAND_4'; 'OUT_GLAC_IMBAL_BAND_5'; 'OUT_GLAC_IMBAL_BAND_6'; 'OUT_GLAC_IMBAL_BAND_7'; 'OUT_GLAC_IMBAL_BAND_8'; 'OUT_GLAC_IMBAL_BAND_9'; 'OUT_GLAC_IMBAL_BAND_10'; 'OUT_GLAC_IMBAL_BAND_11'; 'OUT_GLAC_IMBAL_BAND_12'};
    header = {'YEAR'; 'MONTH'; 'DAY'; 'HOUR'; 'OUT_SWE_BAND_0'; 'OUT_SWE_BAND_1'; 'OUT_SWE_BAND_2';	'OUT_SWE_BAND_3'; 'OUT_SWE_BAND_4'; 'OUT_SWE_BAND_5'; 'OUT_SWE_BAND_6';	'OUT_SWE_BAND_7'; 'OUT_SWE_BAND_8'; 'OUT_SWE_BAND_9'; 'OUT_SWE_BAND_10'; 'OUT_SWE_BAND_11'; 'OUT_SWE_BAND_12'; 'OUT_GLAC_ACCUM_BAND_0'; 'OUT_GLAC_ACCUM_BAND_1'; 'OUT_GLAC_ACCUM_BAND_2'; 'OUT_GLAC_ACCUM_BAND_3'; 'OUT_GLAC_ACCUM_BAND_4'; 'OUT_GLAC_ACCUM_BAND_5'; 'OUT_GLAC_ACCUM_BAND_6'; 'OUT_GLAC_ACCUM_BAND_7'; 'OUT_GLAC_ACCUM_BAND_8'; 'OUT_GLAC_ACCUM_BAND_9'; 'OUT_GLAC_ACCUM_BAND_10'; 'OUT_GLAC_ACCUM_BAND_11'; 'OUT_GLAC_ACCUM_BAND_12'; 'OUT_GLAC_MELT_BAND_0'; 'OUT_GLAC_MELT_BAND_1'; 'OUT_GLAC_MELT_BAND_2'; 'OUT_GLAC_MELT_BAND_3'; 'OUT_GLAC_MELT_BAND_4'; 'OUT_GLAC_MELT_BAND_5'; 'OUT_GLAC_MELT_BAND_6'; 'OUT_GLAC_MELT_BAND_7';	'OUT_GLAC_MELT_BAND_8';	'OUT_GLAC_MELT_BAND_9'; 'OUT_GLAC_MELT_BAND_10'; 'OUT_GLAC_MELT_BAND_11'; 'OUT_GLAC_MELT_BAND_12'; 'OUT_GLAC_OUTFLOW_BAND_0'; 'OUT_GLAC_OUTFLOW_BAND_1'; 'OUT_GLAC_OUTFLOW_BAND_2'; 'OUT_GLAC_OUTFLOW_BAND_3'; 'OUT_GLAC_OUTFLOW_BAND_4'; 'OUT_GLAC_OUTFLOW_BAND_5'; 'OUT_GLAC_OUTFLOW_BAND_6'; 'OUT_GLAC_OUTFLOW_BAND_7'; 'OUT_GLAC_OUTFLOW_BAND_8'; 'OUT_GLAC_OUTFLOW_BAND_9'; 'OUT_GLAC_OUTFLOW_BAND_10'; 'OUT_GLAC_OUTFLOW_BAND_11'; 'OUT_GLAC_OUTFLOW_BAND_12'; 'OUT_GLAC_MBAL_BAND_0'; 'OUT_GLAC_MBAL_BAND_1'; 'OUT_GLAC_MBAL_BAND_2'; 'OUT_GLAC_MBAL_BAND_3'; 'OUT_GLAC_MBAL_BAND_4'; 'OUT_GLAC_MBAL_BAND_5'; 'OUT_GLAC_MBAL_BAND_6'; 'OUT_GLAC_MBAL_BAND_7'; 'OUT_GLAC_MBAL_BAND_8'; 'OUT_GLAC_MBAL_BAND_9'; 'OUT_GLAC_MBAL_BAND_10'; 'OUT_GLAC_MBAL_BAND_11'; 'OUT_GLAC_MBAL_BAND_12'; 'OUT_GLAC_IMBAL_BAND_0'; 'OUT_GLAC_IMBAL_BAND_1';	'OUT_GLAC_IMBAL_BAND_2'; 'OUT_GLAC_IMBAL_BAND_3'; 'OUT_GLAC_IMBAL_BAND_4'; 'OUT_GLAC_IMBAL_BAND_5'; 'OUT_GLAC_IMBAL_BAND_6'; 'OUT_GLAC_IMBAL_BAND_7'; 'OUT_GLAC_IMBAL_BAND_8'; 'OUT_GLAC_IMBAL_BAND_9'; 'OUT_GLAC_IMBAL_BAND_10'; 'OUT_GLAC_IMBAL_BAND_11'; 'OUT_GLAC_IMBAL_BAND_12'};
    columns =length(header) %81 or 82
elseif file_10years(1) == 'c' % cell file is read
%     header = {'YEAR'; 'MONTH'; 'DAY'; 'OUT_PREC'; 'OUT_RAINF'; 'OUT_SNOWF'; 'OUT_EVAP'; 'OUT_EVAP_BARE'; 'OUT_EVAP_CANOP'; 'OUT_TRANSP_VEG'; 'OUT_WDEW'; 'OUT_INFLOW'; 'OUT_RUNOFF'; 'OUT_BASEFLOW'; 'OUT_SWE'; 'OUT_SNOW_DEPTH'; 'OUT_SNOW_CANOPY'; 'OUT_SNOW_COVER'; 'OUT_SNOW_MELT';	 'OUT_SOIL_MOIST_0'; 'OUT_SOIL_MOIST_1'; 'OUT_SOIL_MOIST_2'; 'OUT_GLAC_ACCUM'; 'OUT_GLAC_MELT'; 'OUT_GLAC_SUB'; 'OUT_GLAC_INFLOW'; 'OUT_GLAC_OUTFLOW'; 'OUT_GLAC_WAT_STOR'; 'OUT_GLAC_OUTFLOW_COEF'; 'OUT_GLAC_AREA'; 'OUT_GLAC_MBAL'; 'OUT_GLAC_IMBAL'; 'OUT_GLAC_SURF_TEMP'; 'OUT_GLAC_TSURF_FBFLAG'; 'OUT_GLAC_DELTACC'; 'OUT_GLAC_FLUX'; 'OUT_WATER_ERROR'; 'OUT_AIR_TEMP'; 'OUT_SHORTWAVE'; 'OUT_LONGWAVE'; 'OUT_REL_HUMID'; 'OUT_WIND'; 'OUT_AERO_RESIST'	; 'OUT_AERO_RESIST1'; 'OUT_AERO_RESIST2'};
    header = {'YEAR'; 'MONTH'; 'DAY'; 'HOUR'; 'OUT_PREC'; 'OUT_RAINF'; 'OUT_SNOWF'; 'OUT_EVAP'; 'OUT_EVAP_BARE'; 'OUT_EVAP_CANOP'; 'OUT_TRANSP_VEG'; 'OUT_WDEW'; 'OUT_INFLOW'; 'OUT_RUNOFF'; 'OUT_BASEFLOW'; 'OUT_SWE'; 'OUT_SNOW_DEPTH'; 'OUT_SNOW_CANOPY'; 'OUT_SNOW_COVER'; 'OUT_SNOW_MELT';	 'OUT_SOIL_MOIST_0'; 'OUT_SOIL_MOIST_1'; 'OUT_SOIL_MOIST_2'; 'OUT_GLAC_ACCUM'; 'OUT_GLAC_MELT'; 'OUT_GLAC_SUB'; 'OUT_GLAC_INFLOW'; 'OUT_GLAC_OUTFLOW'; 'OUT_GLAC_WAT_STOR'; 'OUT_GLAC_OUTFLOW_COEF'; 'OUT_GLAC_AREA'; 'OUT_GLAC_MBAL'; 'OUT_GLAC_IMBAL'; 'OUT_GLAC_SURF_TEMP'; 'OUT_GLAC_TSURF_FBFLAG'; 'OUT_GLAC_DELTACC'; 'OUT_GLAC_FLUX'; 'OUT_WATER_ERROR'; 'OUT_AIR_TEMP'; 'OUT_SHORTWAVE'; 'OUT_LONGWAVE'; 'OUT_REL_HUMID'; 'OUT_WIND'; 'OUT_AERO_RESIST'	; 'OUT_AERO_RESIST1'; 'OUT_AERO_RESIST2'};
    columns =length(header) %45 or 46
end


%%%%%
fID = fopen(file_10years, 'r+');
for i=1:6
    line = fgetl(fID);
%     fprintf('%s\n',line);
end
A_10years = fscanf(fID, '%g', [columns Inf])';
fclose(fID);

disp('Number of rows should be:')
rows = size(A_10years,1) %The number of rows in the output file

%%%%%

fID = fopen(file_10yearsWithState, 'r+');
for i=1:6
    line = fgetl(fID);
%     fprintf('%s\n',line);
end
A_10yearsWithState = fscanf(fID, '%g', [columns Inf])';
fclose(fID);

%%%%%

fID = fopen(file_20years, 'r+');
for i=1:6
    line = fgetl(fID);
%     fprintf('%s\n',line);
end
A_20years = fscanf(fID, '%g', [columns Inf])';
fclose(fID);

A_20years = A_20years(end-rows+1:end,:); % cut the lower part of the matrix with appropriate size



%%%%%%%%%%%%%%
% whichColumn = 36; % Define which column you want to compare
nonZeroCounter = 0;

for whichColumn=1:columns

% figure(1)
% plot(A_20years(:,whichColumn), 'b');
% hold on;
% plot(A_10years(:,whichColumn), '--r');
% xlabel('Days')
% ylabel(['Comparision of the ', num2str(whichColumn), 'th column'], 'LineWidth',5)
% title(['20 years vs. 10 years (w/o initial state), (', num2str(whichColumn), ' th column)'])
% legend('20 years','10 years (w/o state)')

if  (sum(abs(A_20years(:,whichColumn))) ~= 0)
    disp(['Normalized abs diff of columns ', num2str(whichColumn), '(', header{whichColumn,1}, ')', ' (20 yrs vs. 10 yrs w/o state)'])
    diff = sum(abs(A_20years(:,whichColumn)-A_10years(:,whichColumn))) / sum(abs(A_20years(:,whichColumn)))
else
    disp(['Zero Column: Abs diff of columns ', num2str(whichColumn), '(', header{whichColumn,1}, ')', ' (20 yrs vs. 10 yrs w/o state)'])
    diff = sum(abs(A_20years(:,whichColumn)-A_10years(:,whichColumn)))
end



% figure(2)
% plot(A_20years(:,whichColumn), 'b');
% hold on;
% plot(A_10yearsWithState(:,whichColumn), '--r');
% xlabel('Days')
% ylabel(['Comparision of the ', num2str(whichColumn), 'th column'], 'LineWidth',5)
% title(['20 years vs. 10 years (w/ initial state), (', num2str(whichColumn), ' th column)'])
% legend('20 years','10 years (w/ state)')

if  (sum(abs(A_20years(:,whichColumn))) ~= 0)
    disp(['Normalized abs diff of columns ', num2str(whichColumn), '(', header{whichColumn,1}, ')', ' (20 yrs vs. 10 yrs w/ state)'])
    diff = sum(abs(A_20years(:,whichColumn)-A_10yearsWithState(:,whichColumn))) / sum(abs(A_20years(:,whichColumn)))
else
    disp(['Zero Column: Abs diff of columns ', num2str(whichColumn), '(', header{whichColumn,1}, ')', ' (20 yrs vs. 10 yrs w/ state)'])
    diff = sum(abs(A_20years(:,whichColumn)-A_10yearsWithState(:,whichColumn)))  
end
disp('---------')


if (diff ~= 0)
    nonZeroCounter = nonZeroCounter + 1;
end


end

nonZeroCounter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fopen()

% load
% fscanf()
% fread()
% fgets()
% fgetl()
% importdata
% sscanf()
% textread()
% textscan()