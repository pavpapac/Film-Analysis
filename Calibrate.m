%%%%%%%%%%%%%%%%%%%%%%% CALIBRATE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Derives the calibration data needed for film
%measurements (single channel) including dose levels, netOD, background, 
%and bad pixel matrix.  INPUT: filename: name of file to save calibration
%data. CHANNEL: 'red', 'green' or 'blue', res: scanning resolution in dpi
%(suggested values: 50 -> 0.5 mm ,127 -> 0.2 mm, 254 - > 0.1 mm)


%% Input parameters 

filename='Ortho120kV_Ref_11000XL_Mar2018.mat';
CHANNEL='green';
res=50;
bkgr=0; % with bkgr (option 1), no bkgr (option 2)
blanc=0; % with blanc scan - option 1 (flood field), no blanc (option 0)

%Provide actual dose levels.
dose=[0 1 2 3 4 5 6 7 8 9 10];
num_of_films=length(dose);

%Define paths and names for scripts and calibration films

addpath(genpath('/Users/pavlos/Dropbox/Personal/MATLAB/FilmQA/'));

filmpath = '/Users/pavlos/Dropbox/Personal/MATLAB/FilmQA/CalibrationFilms_Ortho120kVp/';

filmsets_before = {'Calib120kV_unirr_refl.tif'};

filmsets_after = {'Calib120kV_Arefl.tif'};

film_bkgr={'...'};

film_blank={'Calib120kV_blanc_refl.tif'};

% Select color channel index: 1 (Red), 2 (Green), 3 (Blue)
if strcmp(CHANNEL,'red')
    c=1;
elseif strcmp(CHANNEL,'green')
    c=2;
elseif strcmp(CHANNEL,'blue')
    c=3;
else
    error('Select channel: red, green or blue')
end

%% Convert dpi resolution to mm/pixel
switch(res)
    case 50
        res_mm=0.5;
    case 127
        res_mm=0.2;
    case 254
        res_mm=0.1;
    case 508
        res_mm=0.05;
end

%Background ROI dimensions. Suggested size: 4 cm x 4 cm.
%This corresponds to: 80 (50 dpi), 200 (127 dpi), 400 (254 dpi) and 
% 800 (508 dpi) pixels.

ROI_bkgr=40./res_mm;

%Calibration ROI dimensions. Suggested size: 0.5 cm x 0.5 cm.
%This corresponds to: 10 (50 dpi), 25 (127 dpi), 50 (254 dpi) and 
% 100 (508 dpi) pixels.

ROI_calib=5./res_mm;

% Number of ROIs to be used for calibration. Suggested number: 3-5
num_of_ROIs=3;

%Bad pixel threshold. Suggested values: 0.95-0.99
threshold=0.95;
    
%% Background reading of a scan of a black surface.
%  This will provide us a zero-light intensity value,
%  which should be subtracted by the final film pixel value. We
%  assume that a black surface should always give us pixel value=0.
%  The ROI here could be larger to increase our statistics.
if bkgr==1
    
    scan_bkgr=imread([filmpath film_bkgr{1}]);   
    bkgr=wiener2(scan_bkgr(:,:,c));
    [m_bkgr,s_bkgr] = ROIaverage(bkgr,ROI_bkgr);
else
    m_bkgr=0;
    s_bkgr=0;
end

if blanc==1
    %% Blank scan reading. 
    % This will provide us full-light intensity values for each
    % pixel-detector. We will define a threshold value under which the pixel-detector
    % will be considered ''bad'' and will not be included in the ROI average. Blank
    % scan should be provided as the first film scan before and after
    % irradiation.
    % 

    % Read blank scan and apply Wiener filter.

    scan_blank=imread([filmpath film_blank{1}]);
    blank=wiener2(scan_blank(:,:,c));

    %% Derive a mask of bad pixels.

    figure(2);
    [bad]=findBadPixels(blank,threshold);

    close all;
end

%% Read exposed and unexposed images.

scan_bef=imread([filmpath filmsets_before{1}]);
scan_aft=imread([filmpath filmsets_after{1}]);

%% Apply the wiener filter on regions of 5 x 5 each.

bef=wiener2(scan_bef(:,:,c),[5 5]);
aft=wiener2(scan_aft(:,:,c),[5 5]);

if blanc==1
    %% Apply the bad pixels mask. This will zero all the bad pixels on the images.
    bef=bef.*(1-bad);
    aft=aft.*(1-bad);
end

netOD=zeros(1,num_of_films,'double');
dnetOD=zeros(1,num_of_films,'double');

for j=1:num_of_films

    % Register the exposed to the unexposed film pieces and evaluate ROIs.

    dx=zeros(num_of_ROIs,1,'double');
    dy=zeros(num_of_ROIs,1,'double');

    % Register and choose ROIs for unexposed and exposed films. 

    exposed='no';

    figure(1);

    [bef_reg,x_bef_reg,y_bef_reg] = FilmRegistration(bef);
    [m_bef,s_bef,dx,dy] = ChooseROIs(bef_reg,x_bef_reg,y_bef_reg,dx,dy,num_of_ROIs,ROI_calib,exposed);

    exposed='yes';

    figure(2);

    [aft_reg,x_aft_reg,y_aft_reg] = FilmRegistration(aft);
    [m_aft,s_aft] = ChooseROIs(aft_reg,x_aft_reg,y_aft_reg,dx,dy,num_of_ROIs,ROI_calib,exposed);


    %Calculate netOD and uncertainties (including the background reading) 
    if bkgr==1
        netOD(j)=log10((m_bef-m_bkgr)/(m_aft-m_bkgr));
        netOD(j)=real(netOD(j));          
        dnetOD(j)=(1/log(10)).*sqrt( ((s_bef.^2+s_bkgr.^2) ./ ((m_bef-m_bkgr)^2))+ ((s_aft.^2+s_bkgr.^2)./ ((m_aft-m_bkgr)^2)) );
    else
        netOD(j)=log10((m_bef)/(m_aft));
        netOD(j)=real(netOD(j));          
        dnetOD(j)=(1/log(10)).*sqrt( ((s_bef.^2) ./ (m_bef^2))+ ((s_aft.^2)./ (m_aft^2)) );
    end
    

end

  
%% Reshape arrays in a column format and sort them.

dose=reshape(dose,num_of_films,1);
dose=sort(dose);
netod=reshape(netOD,num_of_films,1);
netod=sort(netod);
dnetod=reshape(dnetOD,num_of_films,1);
dnetod=sort(dnetod);


%% Subtract the netOD of the control film (dose = 0 Gy) and recalculate uncertainty.
dnetod=sqrt(dnetod.^2 + dnetod(1).^2);
netod=netod-netod(1);
 
    
%% Calculate data point weights, calculated as 1/sigma^2. 
%Normalize to the sum of all the weights, so that the magnitude of the
%measured quantity does not dictate the weighting process.

weights=(1./dnetOD).^2;
sum(sum(weights));
weights=weights./sum(sum(weights)); 

%% Now reshape and sort them as well.
weights=reshape(weights,num_of_films,1);
weights=sort(weights,'descend');

close all;  

%% Save calibration data.

save(filename, 'CHANNEL', 'res','dose','netod', 'dnetod', 'weights','m_bkgr','s_bkgr','bad');
disp('Calibration data saved.')

clear variables;
