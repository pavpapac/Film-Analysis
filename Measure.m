%%%%%%%%%%%%%%%%%%%%%% MEASURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUMMARY: function to measure the absolute dose and x/y profiles
% on user-selected or max point. Optionally a dicom image can be imported and
% output and prfile calculations directly compared to measurements.
% Optionally a gamma map and isodose curves are generated on the 2D
% distribution in a user-specified scan range. 

%INPUT: FITTDATA: .mat file containing the fit data
% parameters of the selected protocol (output of script Model.m). CALIBDATA: 
% .mat file containing the calibration data (output of sctipt Calibrate). 
%ROI_dim (mm): ROI dimensions (mm) used for output measurements 
% (side of square sensitive area). For profiles, the resolution is preserved
% (e.g. 127 dpi) on the scanning orientation (x or y), but ROI_dim defines
% the lateral profile width that measurements are averaged on the 
% perpendicular orientation (y or x respectively). The
% dose output is averaged over the (ROI_dim x ROI_dim) square area. 

% POINT. 'origin':  user selects a point on the image by
% clicking on the x and y orientations sequentially that define the origin. 
% It is suggested to mark the x and y orientations on film prior to irradiation. 
% 'max' : the user still selects an (x,y) which is used as an initial point. The code
% then searches the point of maximum dose within the local region (<8 cm x 8 cm)
% surrounding that point.

%POS_ERR = positioning error (mm) due to misplacement of point cursor and/or 
% marker on crosshair. The user may choose to include at this step 
% any other uncertainty (e.g. crosshair to radiation isocenter).  
% The algorithm will choose in total 25 points around the selected pointed
% sampled by a normal distribution with st. dev = pos_err. The reported
% output includes the effect of this positioning error in the 
% uncertainty. NOTE: if POINT = 'max' the pos_err is no longer used. 

%CTRL: Set true if a control film is used. RANGE: The full range (mm) where all
% analysis will occur for both measurements and calculations (-range/2 : range/2)
%Keep this as short as needed to icnrease computational efficiency. 
% Usually 40 mm is sufficient for SRS/ SBRT clinical cases. This will also
% be the range used for the gamma analysis. 

% DICOM. Set tru if a dicom calculation will be provided later for
% comparison with measurements. GAMMA: Set true if a gamma analysis is
% requested between calculated and measured 2D dose distributions. NOTE:
% The GammaCalc script is relatively slow. Thus keep the analysis range as
% short as possible to increase spead. For range=40 mm it should run within
% ~ 1 min. 

%LOCAL: Set to 1 if local gamma is requested or 0 for global gamma. PERC:
%The %dd criterion level, DTA: the distance to agreement level. 

%RES_MM_DCM: The DICOM pixel resolution in mm. Note that the pixel
%resolution does not need to be equal between measurement and calculation
%to perform gamma analysis. 

%INPUT FILES. FILM_BEFORE: The .tif film scan before irradiation, 
%FILM_AFTER: the .tif film scan after irradiation, FILEDICOM: the .dcm
%calculation image. If CTRL==1 then FILM_BEFORE_CTRL and FILM_AFTER_CTRL
%need to be aalso provided. All files need to reside in the same folder or
%in a subfolder of where the Measure script resides. 

% OUTPUT: Absolute dose crossplane (x) and inplane (y) profiles
% of the scanned resolution (e.g. 50,127,254 dpi) with 1 sigma estimated
% uncertainty levels. Also, absolute dose output +-1 sigma reported in
% textbox. if DICOM==1 the dose calculations will be overlaid. 
% A figure with the 25 sampled point will be preserved on screen
% for the user to evaluate the selected points. If GAMMA==1 then a figure
% with the gamma histogram, gamma map and isodose curves will be presented.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%

% Import fit data and calibration data .mat file 

fitdata='fit_green_68hrs_Ref_10000XL_Apr2017.mat';
calibdata='Green_68hrs_Ref_10000XL_Apr2017.mat';

% Set measurement ROI size and positioning error in mm. 

ROI_dim=0.2;
pos_err=0.5;

%Set type of measurement: 'origin' or 'max' and gamma map: 'local' or
%'global'. If no gamma map is needed set 0.  
POINT='max';

% Set true if a control film is used. 
CTRL=1;

% Set the range (mm) that your are interested in performing measurements in
% the films scan. This range will be used later on as the area limit for
% the calculated dose area around the isocenter as [-range/2 range/2]. A
% value of 40 mm or less is suggested for SRS/SBRT beams. 

range=16.2;

%Set DICOM true (1) if you are importing a dicom calculation for comparison and
%provide the filename. Set GAMMA true (1) if you need a gamma map. 
% Choose also if the gamma map should be 'local' (1) or 'global' (0). 

DICOM=0;
GAMMA=0;

% Gamma 'local' (1) or 'global' (0) and gamma criteria levels perc (%) and dta (mm)
local=1;
perc=3; 
dta=3;

%Resolution of the dicom calculations
res_mm_dcm=0.5;


%%%%%%%%%%%%%%%%%% INPUT FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read images before and after irradiation (including control). 
% 
addpath(genpath('/Users/pavlos/Dropbox/Personal/MATLAB/FilmQA/RESEARCH/Films/June2017/'));

% Before irradiation
film_before = 'Glen_TB5_Bef_127dpi.tif';

%After irradiation 
film_after= 'Glen_TB5_127dpi_aft_68hrs.tif';

filedicom='Kulaga-Marzena05mm.dcm';

switch(CTRL)
    case 1
    %Before irradiation (control)
    film_before_ctrl='Glen_TB5_Bef_127dpi.tif';
    %After irradiation (control) 
    film_after_ctrl='Glen_TB5_127dpi_aft_68hrs.tif';
end

%% Load fit and calibration data parameters
fitdata=load(fitdata);
calibdata=load(calibdata);
TYPE=fitdata.TYPE;
channel=fitdata.CHANNEL;
res=fitdata.res;
p=fitdata.p;
q=fitdata.q;
n=fitdata.n;
dp=fitdata.dp;
dq=fitdata.dq;
m_bkgr=calibdata.m_bkgr;
s_bkgr=calibdata.s_bkgr;
bad=calibdata.bad;

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


%% Select color channel index
if strcmp(channel,'red')
    c=1;
elseif strcmp(channel,'green')
    c=2;
elseif strcmp(channel,'blue')
    c=3;
end
    
%%%%%%%%%%%%%%%%%%%%% IMAGE PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Read the films before and after irradiation, select the image 
%channel, apply a wiener filter and correct for bad pixels.

%CONTROL FILM
switch(CTRL)
    case 1
        %Register and select a ROI for unexposed and exposed control films.
        %Suggested ROI size: 1 x 1 cm.
        
        [scan_aft_ctrl,scan_bef_ctrl]=ExtrScans(film_before_ctrl,film_after_ctrl,c,bad);

        num_of_ROIs=1;
        dim=10/res_mm;

        dx=0;
        dy=0; 

        %% Register and rotate films before and after exposure. 

        exposed='no';
        figure(1);
        [bef_reg,x_bef_reg,y_bef_reg] = FilmRegistration(scan_bef_ctrl);
        [m_bef_ctrl,s_bef_ctrl,dx,dy] = ChooseROIs(bef_reg,x_bef_reg,y_bef_reg,dx,dy,num_of_ROIs,dim,exposed);

        exposed='yes';
        figure(1);
        [aft_reg,x_aft_reg,y_aft_reg] = FilmRegistration(scan_aft_ctrl);
        [m_aft_ctrl,s_aft_ctrl] = ChooseROIs(aft_reg,x_aft_reg,y_aft_reg,dx,dy,num_of_ROIs,dim,exposed);

    case 0
        m_bef_ctrl=0;
        s_bef_ctrl=0;
        m_aft_ctrl=0;
        s_aft_ctrl=0;
end

%% Extract and register measurement scans before and after exposure.
[scan_aft,scan_bef]=ExtrScans(film_before,film_after,c,bad);

figure(1);
[scan_aft_reg,x_aft_reg,y_aft_reg] = FilmRegistration(scan_aft);
figure(1);
[scan_bef_reg,x_bef_reg,y_bef_reg] = FilmRegistration(scan_bef);

%% Now extract the measurement area, roi output, crossplane and inplane profiles. 

dx=0;
dy=0;
exposed='yes';
[scan_area_aft,m_roi_aft,s_roi_aft,cross_aft,inpl_aft,dcross_aft,dinpl_aft,pos,dx,dy] ...
    = ExtrMeas(scan_aft_reg,res_mm,x_aft_reg,y_aft_reg,dx,dy,ROI_dim,pos_err,POINT,exposed,range);
exposed='no';
[scan_area_bef,m_roi_bef,s_roi_bef,cross_bef,inpl_bef,dcross_bef,dinpl_bef] ...
    = ExtrMeas(scan_bef_reg,res_mm,x_bef_reg,y_bef_reg,dx,dy,ROI_dim,pos_err,POINT,exposed,range);

%% Now extract dose values from measured values 

[dose2D,dose_roi,dose_cross,dose_inpl,dDtot_roi,dDtot_roi_rel,dDtot_cross,dDtot_inpl] ... 
    = ExtrDose(scan_area_aft,m_bkgr,s_bkgr,m_roi_aft,s_roi_aft,m_aft_ctrl,s_aft_ctrl,cross_aft,inpl_aft,dcross_aft, ... 
    dinpl_aft,scan_area_bef, m_roi_bef,s_roi_bef,m_bef_ctrl,s_bef_ctrl,cross_bef,inpl_bef,dcross_bef,dinpl_bef, ... 
    CTRL,TYPE,p,q,n,dp,dq); 
%% Now extract the dicom 2D dose area, center roi, crossplane and inplane calculations for the
%same plan
switch(DICOM)
    case 1
        [dose2D_dcm,roi_dcm,cross_dcm,inpl_dcm,pos_dcm]=ExtrDICOMcalc(filedicom,res_mm_dcm,range,POINT);
        % Now calculate the gamma map for the two images
        switch(GAMMA)
            case 1
            [gamma,pixels_pass]=GammaMap(dose2D,dose2D_dcm,res_mm,res_mm_dcm,range,local,perc,dta);
        end
       
end

dose_cross=real(dose_cross);
dose_inpl=real(dose_inpl);
dose2D=real(dose2D);
dDtot_cross=real(dDtot_cross);
dDtot_inpl=real(dDtot_inpl);

%%%%%%%%%%%%%%%%%%%%%%% PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2);
subplot(1,2,1);
plot(pos(dose_cross>0),dose_cross(dose_cross>0),'k-');
hold on;
switch(DICOM)
    case 1
        plot(pos_dcm,cross_dcm,'b.');
end
plot(pos(dose_cross>0),dose_cross(dose_cross>0)+dDtot_cross(dose_cross>0),'k--');
plot(pos(dose_cross>0),dose_cross(dose_cross>0)-dDtot_cross(dose_cross>0),'k--');
title('Crossplane','fontsize',15);
xlabel('x (mm)','fontsize',15);
ylabel('Dose (Gy)','fontsize',15);
legend('Film','Calculation',' 1 \sigma');
grid on;

figure(2);
subplot(1,2,2);
plot(pos(dose_inpl>0),dose_inpl(dose_inpl>0),'k-');
hold on;
switch(DICOM)
    case 1
        plot(pos_dcm,inpl_dcm,'b.');
end
plot(pos(dose_inpl>0),dose_inpl(dose_inpl>0)+dDtot_inpl(dose_inpl>0),'k--');
plot(pos(dose_inpl>0),dose_inpl(dose_inpl>0)-dDtot_inpl(dose_inpl>0),'k--');
title('Inplane','fontsize',15);
xlabel('y (mm)','fontsize',15);
ylabel('Dose (Gy)','fontsize',15);
switch(DICOM)
    case 1
    str = {'Measured output:',['(ROI= ' num2str(ROI_dim) ' mm +- ' num2str(pos_err) ' mm) : '], ...
    [num2str(round(dose_roi,2)) '+- ' num2str(round(dDtot_roi,2)) 'Gy ( ' num2str(round(dDtot_roi_rel,2)) ' %)'], ...
    'Calculated output:',[ num2str(round(roi_dcm,2)) ' Gy'], ['% error: ', num2str(round(100.*((dose_roi-roi_dcm)./roi_dcm),2))] };
    case 0
    %str = {'Measured output:',['(ROI= ' num2str(ROI_dim) ' mm +- ' num2str(pos_err) ' mm) : '], ...
    %[num2str(round(dose_roi,2)) '+- ' num2str(round(dDtot_roi,2)) 'Gy ( ' num2str(round(dDtot_roi_rel,2)) ' %)']};
end
%annotation('textbox','Position',[0.3,0.1,0.2,0.2],'String',str,'fitBoxToText','on');
grid on;

if (DICOM==1 && GAMMA==1)
    figure(3);
    subplot(1,3,1);
    histogram(gamma);
    title([num2str(round(pixels_pass,2)) ' % pixels pass the gamma test'],'fontsize',15);
    xlabel('gamma','fontsize',15)
    ylabel('# pixels','fontsize',15)
    
    figure(3);
    subplot(1,3,2);
    iso=find(pos==0);
    dose_iso=dose2D(iso,iso);
    isolevels=dose_iso.*[0.1 0.2 0.3 0.5 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5];
    contour(pos,pos,dose2D,isolevels,'-');
    hold on;
    iso_dcm=find(pos_dcm==0);
    dose_iso_dcm=dose2D_dcm(iso_dcm,iso_dcm);
    isolevels_dcm=dose_iso_dcm.*[0.1 0.2 0.3 0.5 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5];
    contour(pos_dcm,pos_dcm,dose2D_dcm,isolevels_dcm,'--');
    grid on;
    colormap default;
    colorbar;
    legend('Measured', 'Calculated');
    title('Isodose map','fontsize',15)
    xlabel('x (mm)','fontsize',15)
    ylabel('y (mm)','fontsize',15)
    
    figure(3);
    subplot(1,3,3);
    imagesc(pos_dcm,pos_dcm,gamma);
    colormap default;
    colorbar;
    title('Gamma map','fontsize',15)
    xlabel('x (mm)','fontsize',15)
    ylabel('y (mm)','fontsize',15)
    
end
   
cross=dose_cross./max(dose_cross);
cross=(flip(cross)+cross)./2;

clearvars -except x dose_cross cross dDtot_cross crossA crossB crossC crossD crossE dose_crossA dDtot_crossA dose_crossB dDtot_crossB dose_crossC dDtot_crossC dose_crossD dDtot_crossD dose_crossE dDtot_crossE;
%disp(str) 
