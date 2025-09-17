function [SC,T] = combineSturmJohnstonSnowClasses(path_Sturm,path_Johnston,path_elevation)
%combineSturmJohnstonSnowClasses Produces a combined snow classification
%system using the Sturm and Liston 2021 snow classes from 'Revisiting the
%Global Seasonal Snow Classification' and that from Johnston et al., 2023
%'Global Snow Seasonality Regimes...' using ~1km resolution data.
%
% RELEVANT LITERATURE AND DATASETS (linked below):
%
%   1. Sturm, M. & Liston, G. E. Revisiting the Global Seasonal Snow 
%      Classification: An Updated Dataset for Earth System Applications. 
%      Journal of Hydrometeorology 22, 2917–2938 (2021).
%
%   2. Johnston, J., Jacobs, J. M. & Cho, E. Global Snow Seasonality 
%      Regimes from Satellite Records of Snow Cover. Journal of 
%      Hydrometeorology 25, 65–88 (2023).
%
% INPUTS:
%
%   path_Sturm - file path to Sturm and Liston global (GL) 30 arc second
%   dataset. Available for download at NSIDC: 
%   < https://nsidc.org/data/NSIDC-0768/versions/1 >
%
%   path_Johnston - file path to Johnston et al., 2023 snow seasonality
%   classification (SSC) dataset (0.01 degree grid). Available for download
%   at NSIDC: < https://nsidc.org/data/nsidc-0791/versions/1 >
%
%   path_Elevation - file path to GMTED2010 30 arc second elevation data.
%   Available for download: < http://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/Grid_ZipFiles/mn30_grd.zip >
%
% OUTPUTS:
%
%   SC - the snow classification as a 2D grid with 30 arc second
%   resolution. This is also saved as a geotiff file at the same path as
%   the input Sturm and Liston dataset (path_Sturm)
%   
%   T - a table legend with the numeric value and its corresponding snow
%   class. This table also includes recommended color codes for
%   visualization.
%
% METHODS:
%
% To produced the combined snow classification, this script executes the
% following steps.
%
%   1) Resamples the snow seasonality averages classes from the Johnston et
%   al., 2023 dataset (1/100 degree grid) to match the Sturm and Liston
%   2021 seasonal snow classification dataset (1/120 degree grid) using
%   nearest neighbor interpolation.
%
%   2) Renumbers the Sturm and Liston classification to values
%   10 (ephemeral), 20 (prairie), 30 (maritime), 40 (montane forest), 
%   50 (boreal forest), 60 (tundra), and 70 (snow/ice) then adds the
%   average snow seasonality classification values (0-4). Where 0 = regions
%   with no observed snow cover in all snow years and 4 = permanent snow
%   cover observed in all snow years. Johnston et al., 2023 classes
%   represent no snow (0), ephemeral (1), transitional (2), seasonal (3),
%   and perennial (4) snow. The Sturm and Liston 2021 dataset was 
%   re-numbered to better reflect seasonality, in which lower values tend 
%   to be in warmer or ephemeral snow regions (10,20,30), and the higher 
%   values (50,60,70) reflect very cold climates.
%
%   3) Classes are rounded to the nearest integer value, thus must take on
%   1 of 35 discrete values (10,11,12,13,14,20,21,...70,71,72,73,74).
%
%   4) Classes known to not occur in nature (e.g., 70: no snow snow/ice)
%   and those will very small spatial exents are combined with the nearest
%   logical class. Resulting in 18 snow classes. The combination methods 
%   are detailed in the following script and in the accompanying manuscrupt
%   (DOI: TBD).
%
%   5) Perform artifact filtering using the following steps:
%       a. identify all unique connected areas which contain ephemeral
%          ephemeral (11) or transitional ephemeral (12) snow classes
%       b. iterate through each connected group and calculate the mean
%          elevation of the interior area within the class, as well as the
%          mean elevation and the snow classes of the pixels directly 
%          bordering the connected areas
%       c. if the mean elevation of the surrounding boundary (+10m) is
%          greater than that of the interior area and >90% of the pixels in
%          the surrounding area are defined as no snow (10) or water (0),
%          the interior region is set as no snow (10). This is intended to
%          identify regions of similar or lesser elevations that have snow
%          classes of higher seasonality than the surrounding areas.
%          Primarily, this identifies and removes the effect of salt flats 
%          and lake beds on the classification. Which may result in
%          inaccurate representations of the snow seasonality due to 
%          similar Normalized Difference Snow Index (NDSI) of these areas
%          to snow cover.
%   
%   6) Step 5 is repeated as above, but instead used to identify 
%   transitional ephemeral (12) areas that are surrounded by ephemeral (11)
%   and water (0) classes. This removed additional water body influences in
%   ephemeral snow areas.
%
%   7) Produce a legend .csv table relating the numeric values to each of
%   the 18 snow classes. Also, includes recommended color values for
%   visualizing the the snow class data, which is saved as a geotiff file 
%   to the same path as the input 'path_Sturm' data.
%
% GENERAL NOTES:
%
% This script is designed to work with datasets of the noted resolutions
% above, it is not adaptive to inputting snow classification datasets for
% different spatial resolutions or with varied domain coverage. If using
% this script to produce datasets at other resolutions, please modify the
% code accordingly to suit the input snow class and elevation datasets.
%
% Testing on done on a Mac M1 chip with 16 GB of RAM using MATLAB R2024b.
% Total processing time was approximately 2 minutes.

%time script execution
tt = tic;

%load in Sturm seasonal snow classification data (1/120 degree grid)
disp('Loading in snow classification data...');
[SC1,R1] = readgeoraster(path_Sturm);

%load in the climatology average snow class, considering all years in the
%published Johnston et al., 2023 dataset
SC2 = ncread(path_Johnston,'snow_class_climatology');
SC2 = flipud(rot90(SC2)); %reorient so north is up and east is right

%load in 1/120 degree (30 arc second) elevation data (only for 90S to 84N)
[E,~] = readgeoraster(path_elevation); %does not have georef info (other finer grids do)

%set all elevation values as NaN above 84N to match grid extent to that of
%the Sturm and Liston 1/120 degree classification
Em = NaN(size(SC1));
Em((size(SC1,1) - size(E,1)) + 1:end,:) = E;

%Johnston mapping covers -180 to 180 longitude and -90 to 90 latitude,
%latitude and longitude values are included as cell center coordinates of
%each 0.01 degree grid cell, thus they max at 0.005 short of these values
%(e.g., max lat = 89.995). The Sturm global classification also covers the
%same geographic extent.
% Thus, we simply resample the nearest data point in the 1/100 degree 
% grid used by Johnston using the 1/120 degree grid (30 arc second) used by
% Sturm
SC2 = imresize(SC2,size(SC1),"nearest");

%we renumber below when merging classes, ephemeral = 10, prairie = 20,
%maritime = 30, montane forest = 40, boreal forest = 50, tundra = 60,
%snow/ice = 70
SC1(SC1 == 4) = 10;
SC1(SC1 == 5) = 20;
SC1(SC1 == 3) = 30;
SC1(SC1 == 6) = 40;
SC1(SC1 == 2) = 50;
SC1(SC1 == 1) = 60;
SC1(SC1 == 7) = 70;

%update time elapsed
toc(tt);

%create combined classification by adding the decimal average snow class
%values from Johnston (SC2, 0-4) to the above Sturm values (SC1, 10-70)
disp('Combining classifications...');
SC = single(SC1) + single(SC2);

%set water/non-classified areas in either product to NaN
SC(SC1 == 8 | isnan(SC2)) = NaN;

%round to nearest integer, for final classification, NaNs become 0, which
%is considered the fill/no data value
SC = uint8(round(SC));

%combine 35 classes based on extent and similarity, used to reduce from 35
%classes down to 18 (plus 0 == water).
% All classes defined as perennial ice and snow in the Sturm classification
% (70) were retained as perennial ice and snow, regardless of the snow 
% cover classification value. The proposed combination defaults to the 
% Sturm and Liston 2021 classification in most cases with the assumption 
% that there are imperfections in MODIS-based SCA estimates due to cloud 
% cover, vegetation effects, and surface features (e.g., sand flats) that 
% may influence snow cover mapping accuracy. Classification combination
% rules are detailed as follows and were informed by an analysis into the
% spatial coverages of each of the 35 snow classes:
SC(SC >= 70) = 74; %70,71,72,73,& 74 = 74: all perennial ice and snow (70,71,72,73,74)
SC(SC == 64) = 63; %63,64 = 64: perennial tundra (64) small extent was merged with seasonal tundra, large extent (63)
%62 = 62: no change to transitional tundra (62) extent
SC(SC == 60) = 61; %60,61 = 61: no snow tundra (60) very small extent, merged with ephemeral tundra (61)
SC(SC == 54) = 53; %53,54 = 53: perennial boreal forest (54) very small extent (~300 km2) was merged with seasonal boreal forest (53)
SC(SC == 50 | SC == 51) = 52; %50,51,52 = 52: no snow (50) or ephemeral boreal forest (51) had very small extents and do not realistically exist in nature, combined with transitional boreal forest (52)
SC(SC == 44) = 43; %43,44 = 43: perennial montane forest (44) covered only 27 km2, merged with seasonal montane forest (43) which covered >500,000 km2
%42 = 42: no change to transitional montane forest (42) extent
SC(SC == 40) = 41; %40,41 = 41: no snow montane forest (40) covered a very small extent, was merged with ephemeral montane forest (41)
SC(SC == 34) = 33; %33,34 = 33: perennial maritime (34) covered a very small exent, merged with seasonal maritime (33)
%32 = 32: no change to transitional maritime (32) extent
SC(SC == 30) = 31; %30,31 = 31: no snow maritime (30) was over 50,000 km2 in extent, but was merged with the Sturm ephemeral maritime (31) class to ensure consistency in class boundaries, low snow cover in these areas may be due to dense tree cover
SC(SC == 24) = 23; %23,24 = 24: perennial prairie (24) <800 km2, merged with seasonal prairie (23)
%22 = 22: no change to transitional prairie (22) extent
SC(SC == 20) = 21; %no snow prairie (20) had substantial spatial coverage ~300,000 km2, but merged with ephemeral prairie (21) to maintain Sturm class boundaries
%special cases:
% there was a very small extent of cases in which Sturm ephemeral 
% overlapped with seasonal (13) or perennial (14) areas as classified in 
% the Johnston et al., 2023 dataset (<7000 km2). These were noted to be 
% primarily salt flats and areas around lake boundaries and are set as 
% ephemeral-ephemeral (11). Filtering steps below were applied to change
% false transitional-ephemeral classes to ephemeral, then subsequently
% false ephemeral-ephemeral classes to no snow using the filtering approach
% detailed below.
SC(SC == 13 | SC == 14) = 11;
%12 = 12: no change to transitional ephemeral (12) extent
%10 = 10: no change to no snow ephemeral (10) extent

%update time elapsed
toc(tt);

%apply final filtering step to remove lake/salt flat effects which cause
%false snow classifications due to similar reflective properties to snow
disp('Filtering classes...');

%(1) Filter out false ephemeral areas surrounded by no snow areas
%identify all ephemeral areas
ephemeral_idx = (SC == 11 | SC == 12);

%identify all connected ephemeral areas and the number of pixels within them
BW = bwlabel(ephemeral_idx);
dims = size(BW);

%get linear indices of all unique connected areas
CC = bwconncomp(BW);

%loop through each connected area
for i = 1:(CC.NumObjects)

    %get linear indices of area
    I1 = CC.PixelIdxList{i};

    %convert linear indices to row,col
    [row,col] = ind2sub(dims,I1);

    %get linear indices of all neighboring pixels
    I2 = getLinearIndicesNeighbors(dims,row,col);

    %remove the indices from the original input
    I2 = setdiff(I2,I1);

    %get snow classes and elevation of the neighboring region
    SCdil = SC(I2);
    Edil = single(Em(I2));

    %get elevation of interior regions
    Eint = single(Em(I1));

    %if the mean elevation of the dilated boundary +10m is above that of the
    %interior ephemeral area, and the surrounding classifications are
    %only >90% water (0) and no snow (10), set as no snow
    if mean(Edil)+10 > mean(Eint)
        if sum(SCdil == 10 | SCdil == 0) > length(SCdil)*0.9
            SC(I1) = 10;
        end
    end
end

%(2) Filter out transitional ephemeral areas surrounded by ephemeral areas
%of the same approximate elevation
%identify all ephemeral areas
trans_ephemeral_idx = (SC == 12);

%identify all connected transitional ephemeral areas and the number of pixels within them
BW = bwlabel(trans_ephemeral_idx);
dims = size(BW);

%get linear indices of all unique connected areas
CC = bwconncomp(BW);

%loop through each connected area
for i = 1:(CC.NumObjects)

    %get linear indices of area
    I1 = CC.PixelIdxList{i};

    %convert linear indices to row,col
    [row,col] = ind2sub(dims,I1);

    %get linear indices of all neighboring pixels
    I2 = getLinearIndicesNeighbors(dims,row,col);

    %remove the indices from the original input
    I2 = setdiff(I2,I1);

    %get snow classes and elevation of the neighboring region
    SCdil = SC(I2);
    Edil = single(Em(I2));

    %get elevation of interior regions
    Eint = single(Em(I1));

    %if the mean elevation of the dilated boundary +10m is above that of the
    %interior transitional ephemeral area, and the surrounding classifications are
    %>90% water (0) and ephemeral (11), set as ephemeral
    if mean(Edil)+10 > mean(Eint)
        if sum(SCdil == 11 | SCdil == 0) > length(SCdil)*0.9
            SC(I1) = 11;
        end
    end
end

%update time elapsed
toc(tt);

%update status and create savepath
idx = strfind(path_Sturm,'/');
savepath = path_Sturm(1:idx(end));
disp(['Saving output files to ' savepath]);

%produce table including a legend and hexidecimal colorcode for each class
U = unique(SC);
T = table();
T.snow_class_ID = U;
T.snow_class = cell(size(U));
T.HTML_color = cell(size(U));
T.RED = NaN(size(U));
T.GREEN = NaN(size(U));
T.BLUE = NaN(size(U));
snow_classes = {'no data/water','no snow','ephemeral ephemeral','transitional ephemeral',...
    'ephemeral prairie','transitional prairie','seasonal prairie',...
    'ephemeral maritime','transitional maritime','seasonal maritime',...
    'ephemeral montane forest','transitional montane forest','seasonal montane forest',...
    'transitional boreal forest','seasonal boreal forest',...
    'ephemeral tundra','transitional tundra','seasonal tundra','perennial snow and ice'};
HTML_colors = {'#ffffff','#fff5d2','#e6e673','#b4b400','#f0be91','#f07828',...
    '#be5000','#e68c8c','#dc2828','#960000','#b4d7aa','#6ec85f','#289600',...
    '#648cc8','#0050c8','#dcd2e6','#aa73e6','#7800e6','#9b9b9b'};
RED = [255,255,230,180,240,240,190,230,220,150,180,110,40,100,0,220,170,120,155];
GREEN = [255,245,239,180,190,120,80,140,40,0,215,200,150,140,80,210,115,0,155];
BLUE = [255,210,115,0,145,40,0,140,40,0,170,95,0,200,200,230,230,230,155];
L = length(U);
for row = 1:L
    T.snow_class{row} = snow_classes{row};
    T.HTML_color{row} = HTML_colors{row};
    T.RED(row) = RED(row);
    T.GREEN(row) = GREEN(row);
    T.BLUE(row) = BLUE(row);
end

%save the output as a geotiff and table to the same path as the original sturm
%classification
writetable(T,[savepath '18snowclasses_Sturm2021_Johnston2023_merged.csv']);
geotiffwrite([savepath '18snowclasses_Sturm2021_Johnston2023_merged_01km_30.0arcsec.tif'],SC,R1);

%displau processing time
disp('Processing completed.');
toc(tt);

end