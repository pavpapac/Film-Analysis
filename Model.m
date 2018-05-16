%%%%%%%%%%%%%%%%%%%%% MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Fits a parametric model to the data obtained with
%Calibrate.m script. INPUT: calibdat: .mat file with calibration data, 
%TYPE: 'rational' - > (p*x)/(1+q*x^n) or 'polynomial' -> p*x+q*x^n function
% n: problem-dependent parameter tuned by the user based on uncertanity vs 
% accuracy analysis, filename: name of file to save fit data. 

%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%

calibdata='Green_20hrs_Ref_10000XL_Apr2017.mat';
TYPE='polynomial';
n=2.6; % usually about 2.5-2.7 for reflection @ green 
filename='fit_Green_20hrs_Ref_10000XL_Apr2017.mat';

%% LOAD DATA

data=load(calibdata);
dose=data.dose(1:12);
netod=data.netod(1:12);
dnetod=data.dnetod(1:12);
weights=data.weights(1:12);
CHANNEL=data.CHANNEL;
res=data.res;



%% SELECT MODEL: polynomial or rational. Power parameter remains open for user to tune. 

if strcmp(TYPE,'rational')
    %typical value: n=0.5
    f=fittype('(p*x)/(1+q*x^n)','problem','n');
    [c,gof]=fit(netod,dose,f,'Weight',weights,'problem',n);
elseif strcmp(TYPE,'polynomial')
    %typical value: n=2.5
    f=fittype('p*x+q*x^n','problem','n');
    [c,gof]=fit(netod,dose,f,'Weight',weights,'problem',n);
else
    error('Choose one of the following: rational or polynomial');
end

% Extract fitting uncertainties.

[p,q,dp,dq]=FitParametersErrors(c,gof);
n=c.n;

num_dose=length(dose); 
dose=dose(1:num_dose);
netod=netod(1:num_dose);
dnetod=dnetod(1:num_dose);
x=netod;
dx=dnetod;

if strcmp(TYPE,'rational') % Dfit=p*x/(1+q*x^n)
    
    %First calculate the dose estimated by the fitting function
    
    dose_fit=p*x./(1+q*x.^n);
    
    %Now calculate the error between fit and doe delivered.
    
    err_dose=100.*(dose_fit-dose)./dose;
    
    %Now calculate the partial derivatives
    
    dD_dp=x./(1+q*x.^n);
    dD_dq=(p*(x.^(n+1)))./((1+q*(x.^n)).^2);
    dD_dx=((p*(1+q*x.^n)-n*p*q*x.^n))./((1+q*x.^n).^2);

    
    
elseif strcmp(TYPE,'polynomial') %Dfit=p*x+q*x.^n
    
    %First calculate the dose estimated by the fitting function
    
    dose_fit=p*x+q*x.^n;
    
    %Now calculate the absolute relative difference with dose delivered
    
    err_dose=100.*(dose_fit-dose)./dose;
    
    %Now calculate the partial derivatives
    
    
    dD_dp=x;
    dD_dq=x.^n;
    dD_dx=p+n*q*x.^(n-1);
    
    
end

%Now calculate the fit error
    
dDfit=sqrt( (dD_dp.*dp).^2 + (dD_dq.*dq).^2);

%Now calculate the experimental error
dDexp=sqrt( (dD_dx.*dx).^2 );
    
%Now calculate the total error
dDtot=sqrt( (dD_dp.*dp).^2 + (dD_dq.*dq).^2 + (dD_dx.*dx).^2); 

%Now normalize to the dose level
    
dDfit=100*(dDfit./dose);
dDexp=100*(dDexp./dose);
dDtot=100*(dDtot./dose);


% Now calculate the number of dose points that agree within the estimated 
% 1 sigma level. Assuming the data points follow a normal distribution we
% should expect an accurate model to agree ~ 68.3 % of the times with the
% estimated uncertainty. 

points_1sigma=abs(err_dose)<dDtot;
num_1_sigma=uint8(100.*sum(points_1sigma)./length(dose));
num_1_sigma=int2str(num_1_sigma);

%% PLOTS

%First plot the calibration dose points and the fitted curve. 

figure
subplot(2,2,[3,4]);
plot(netod,dose,'k*');
hold on;
plot(c,'k-');
legend('Data points', 'fitted curve','Location','best');
title('Film calibration curve','fontsize',15);
xlabel('netOD','fontsize',15);
ylabel('Dose (Gy)','fontsize',15);
%ylim([min(dose) max(dose)]);
%xlim([min(netod) max(netod)]);

%Then, plot the total, experimental and fitting uncertainties
subplot(2,2,1);
plot(dose,dDtot,'k^-');
hold on;
plot(dose,dDexp,'rs--');
plot(dose,dDfit,'bo-');
legend('Total', 'Exp', 'Fit', 'Location', 'best');
title('Uncertainty analysis', 'FontSize',15,'LineWidth',2);
xlabel('Dose (Gy)','FontSize',15);
ylabel('\deltaD/D (%)','FontSize',15);
xlim([min(dose) max(dose)+0.1*max(dose)]);
grid on;

%Now plot the relative absolute error as a function of dose and using as errorbars
%the calculated precision for that dose level. This curve serves as an accuracy
%check of the model relative to our measurements. 

subplot(2,2,2)
str = [num_1_sigma,' % of data points < 1 \sigma '];
errorbar(dose,err_dose,dDtot,'ko-');
legend(str,'location','best');
title('Accuracy analysis (residuals)' , 'FontSize',15,'LineWidth',2);
xlabel('Dose (Gy)','FontSize',15);
ylabel('(Dfit-D)/D (%)','FontSize',15);
xlim([min(dose) max(dose)+0.1*max(dose)]);
grid on;

 %% Save data
save(filename, 'TYPE','CHANNEL','res','p','q','n','dp','dq');
disp('Fit data saved.');

clear all;