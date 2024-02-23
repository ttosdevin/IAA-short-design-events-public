%% Short design events for waves and wind (MLER and CMLERs (CRRWs))
% Script to produce combined response conditioned wave and wind events for FOWTs.
% The inputs are RAOs calculated from OpenFAST with a wave only white noise simulation
% and a wind only simulation (U= 50m/s, TI = 11%) with the turbine parked.
% The mean surge offset in the sea state is also needed, a step wind
% simulation could be used to estimate this or a frequency domain model such as RAFT. 

% The MLER is a single event, the constrained MLERs (CMLERs) can be used to
% predict the characteristic value by taking the mean of a sample of at
% least 6 cases. 

% This script is a bit different (and improved) from the one used to produce the cases from
% the report 'Impact accelerator award - Short design events' as it uses an analytical approach to get the combined
% distribution. This requires the use of the Gausian distribution of the process rather than the rayleigh distributed peaks but the results are very similar.

% Contact: tom.tosdevin@postgrad.plymouth.ac.uk 

%% Script 

% 1. Define cases and load wave RAOs
% 2. Load and reprocess wave RAOs, create response spectrum, estimate target response EVD due to waves
% 3. Load and reprocess wind RAOs, create response spectrum, estimate target response EVD due to wind
% 4. Estimate the relative importance of the waves relative to wind and calculate a weighting
% 5. Estimate the target combined response extreme value distribution from waves and wind  
% 6. Produce the wave cases
% 7. Produce the wind cases
% 8. Save data


%% 1. Define cases
clear
close all
clc
tic;

numsimMC = 0; % number of montecarlo simulations to check the target EVDs
MCTS = 0.2; % Time step used in the MC simulations (unused if numsimMC = 0)
EX = 1; % EX = 1 is extreme sea state (U = 50 m/s), 3 is U = 29 m/s
Responsetype = 1; % 1 = Tower base bending moment, 2 = FT1 (mooring load), 3 = Nacelle aceleration in x, 4 = Pitch
Device = 1; % 1 = Windcrete
saveq = 1; % 1 = save, 2 = dont save
plots = 2; % 1 = plot, 2 = dont plot
RIC = 100; % value below which cases are split, if always want to split, use 100.
CRRWnum = 50; % number of CMLERs to generate for average used in place of MLER for wind. 
CRRWnumSave = 10; % number of CMLERs to save
ddt = 0.025; % time step used in output time sereis, should be the same as will be used in OpenFAST.
hourexposureLFMLER = 1/6; % Exposure time of wind in hours, 10 minute repeating used in report
hourexposureWF = 1; % Exposure time of waves (hrs)
ttMLER = -400:ddt:100; % time series time steps with the extreme response occuring at t = 0s.
scalefactor = 1; % scale of the model, 1 if full scale
 d = 200; % water depth

% create name
if EX == 1
    pref1 = 'EX_U50';
end
if EX == 3 
    pref1 = 'mod_U29';
end
if Device == 1
    pref2 = 'WC';
    pref4 = '_08'; % target response amplitude percentile
end

if Responsetype == 1
    pref3 = 'TwrMy';
end
if Responsetype == 2
    pref3 = 'FT1';
end
if Responsetype == 3
    pref3 = 'Nxa';
end
if Responsetype == 4
    pref3 = 'P';
end

fnamMLERWind = strcat('MLER_',pref1,'_wind_',pref2,'New_',pref3,pref4); % wind
fnamMLER = strcat('MLER_',pref1,'_wave_',pref2,'New_',pref3,pref4); % 
fnamCMLERWind = strcat('CMLER_',pref1,'_wind_',pref2,'New_',pref3,pref4); % wind
fnamCMLER = strcat('CMLER_',pref1,'_wave_',pref2,'New_',pref3,pref4); % 

fnamCMLERWindmean = strcat('Mean_CMLER_',pref1,'_wind_',pref2,'_',pref3,pref4); % wind
fnamEVD = strcat('EVD_',pref1,'_wave_',pref2,'_',pref3,pref4); % 

if Device == 1 
percentile = 0.8; % target percentile, recomend 0.8 for spars and TLPs, 0.99 for semi-subs and barges
end
percentileComp099 = 0.99;
percentileComp05 = 0.5;

% mean surge offset due to constant wind, need to write generic interpolated
% script for this based on step wind run. 
if Device == 1
        if EX == 1
XX = 5; % U = 50m/s % this needs predicting from a step wind simulation
        else
XX = 1.7; % U = 29m/s % this needs predicting from a step wind simulation
        end
end

if Responsetype == 2
    XX = XX*1.3; % inflate the target focus location for the wave for the mooring load case. This is because a large surge position usually leads to the extreme mooring response.
end

% Load the RAO wave data. 
load('WhiteNoise_10000_RAOs')

%% 2. Load and reprocess wave RAOs, create response spectrum, estimate target response EVD due to waves
 
% create wave spectrum
if EX == 3
TP = 9;
HS = 5.1;
U = 29 ;
end
if EX == 1
TP = 14.2;
HS = 10.7;
U = 50 ;
end
gamma = 3.3;
  
%Name spectral type to compare to
specType = 'JONSWAP_1';

%Define spectral parameters of run
param = [HS TP gamma];
% orig used
val1 = 1/(TP*2); % fmin
val2 = 1; % fmax
RT1 = hourexposureWF*60*60; % repeat time 

% define angular frequencies
fHz = [ceil(val1*RT1):1:floor(val2*RT1)]/RT1; % frequencies to cover
f = fHz*(2*pi); %convert to angular frequency
wvec = f;
Period = (2*pi)./wvec;

fc = linspace(f(2)-f(1),f(2)-f(1),length(f)); % delta f

[Sf] = waveSpectrum(fHz, param, specType);
SWave = Sf/(2*pi);
fname = [];

% set these with reasonable values of I to avoid large file size plots etc.
% I should cover a range coving the extreme response
if Responsetype == 4
RAO = RAOP;
RAOpha = PhasePitch;
I = 0:0.01:30; % response values (deg)
end
if Responsetype == 1
RAO = RAOTMy;                  
RAOpha = PhaseTMy;
I = 1000:1000:10000000; % response values (Nm)
end
if Responsetype == 2
RAO = RAOLC1;               
RAOpha = PhaseLoadcell1;
I = 10000:1000:100000000; % response values (N)
end
if Responsetype == 3
RAO = RAONx;                  
RAOpha = PhaseNx;
 I = 0:0.01:20; % response values (m/s^2)
end

SHz = SWave*(2*pi);

%% Resample so that spectrum and RAO frequencies match up

%resample RAO at frequencies of S
fRAO = 1./fRAO; % as fRAO in Hz
fRAO(isinf(fRAO)|isnan(fRAO)) = fRAO(2)+(0.01*fRAO(2));
% fRAO = 1./fRAO; % as fRAO in Hz
tsin = timeseries(RAO,fRAO);
tsout = resample(tsin,Period);
RAOsampledWave = flip(squeeze(tsout.data));

tsin = timeseries(RAOpha,fRAO);
tsout = resample(tsin,Period);
phaseofinterestWave = flip(squeeze(tsout.data));
SrWave = (SWave).*(RAOsampledWave.^2); % convert to response spectrum. 
fWave = wvec;

% plots of RAOs vs sampled
% figure
% hold on
% plot(fRAO,RAO)
% plot(Period,(RAOsampledWave))
%  
% figure
% plot(Period,SWave)
% hold on
% yyaxis right
% plot(Period,(SrWave))
%% Calculate spectral moments and amplitudes


% S in frequency components (for M0)
b = (SrWave').*fc;
% for M1
b1 = (SrWave').*fc.*(f);
%M2
b2 = (SrWave').*fc.*(f.^2);
% M4   
b4 = (SrWave').*fc.*(f.^4);
WFbj = sqrt(2*SrWave'.*fc);

% spectral moments 
M0Wave = sum(b,2);
M1Wave = sum(b1,2);
M2Wave = sum(b2,2);
M4Wave = sum(b4,2);
MM0Wave = sqrt(M0Wave);
b = [];
b1 = [];
b2 =[];
b4 = [];
% calculate target amplitude
% from moments of response spectrum
ET = (((hourexposureWF*60)*60)/(sqrt(scalefactor))); % Exposure time
ee = sqrt(1-(((M2Wave)^2)/(M0Wave*M4Wave))); %bandwidth parameter
n= (1/(4*pi))*((1+sqrt(1-(ee^2)))/sqrt(1-(ee^2)))*sqrt(M2Wave/M0Wave)*(ET); % (2.12) expected number of response peaks in exposure time
TEV = sqrt(2*log((2*n*sqrt(1-(ee^2)))/(1+sqrt(1-(ee^2))))); % Target extreme value, 

nWF_old = n;

% predict wave response EVD
FFmax =[];
eeee = [];
sig = 1;
muWF = TEV*sig;
m = round(1/(normcdf(0,muWF,sig)));

% Rayleigh peaks
% for i = 1:length(I)
%     eeee(i) = I(i)/MM0Wave;
% FFmax(i) = n*((I(i)/((MM0Wave^2)))*exp(-(I(i)^2)/(2*(MM0Wave^2)))*(1-exp(-(I(i)^2)/(2*(MM0Wave^2))))^(n-1));
% end

% Gausian process
% for i = 1:length(I)
%     eeee(i) = I(i)/MM0Wave;
% FFmax(i) = m*((1/((MM0Wave*sqrt(2*pi))))*exp(-(I(i)^2)/(2*(MM0Wave^2))))*(normcdf(eeee(i),0))^(m-1);
% end
        eeee = I/MM0Wave;
        FFmax = m*((1/((MM0Wave*sqrt(2*pi)))).*exp(-(I.^2)/(2*(MM0Wave^2)))).*(normcdf(eeee,0)).^(m-1);

tarWF = cumtrapz(eeee*MM0Wave,FFmax); % get CDF from PDF (extremes)
if plots == 1
plot(eeee*MM0Wave,tarWF,'b-.') 
hold on
end

tarpksWF = nthroot(tarWF,n); % CDF of peaks rather than extremes by taking nth root
WFXaxis = eeee*MM0Wave;

[minValue,closestIndex] = min(abs(tarWF-percentileComp099'));
MWF099 = eeee(closestIndex)*MM0Wave; % Target amplitude of response at 99th percentile
[minValue,closestIndex] = min(abs(tarWF-percentileComp05'));
MWF05 = eeee(closestIndex)*MM0Wave; % Target amplitude of response at 50th percentile
[minValue,closestIndex] = min(abs(tarWF-percentile'));
MWF = eeee(closestIndex)*MM0Wave; % Target amplitude of response


%% 3. Load and reprocess wind RAOs, create response spectrum, estimate target response EVD due to wind


% Load wind RAOs

load('RAO_U50_10000_RAOs.mat')

% define wind spectrum
TI = 0.11; % turbulence intensity
sig = TI*U;
RT1 = 10*60; % repeat time = 10 mins
RT = 10*60; % repeat time 

val1 = 1/RT1; %  fmin 
val2 = 0.3; % fmax

% define angular frequencies
fHz = [ceil(val1*RT1):1:floor(val2*RT1)]/RT1; % frequencies to cover (Hz)

%%
f = (2*pi).*fHz; % angular freq
fc = linspace(f(2)-f(1),f(2)-f(1),length(f)); % delta f

Sf = ((sig^2)*(4*340.2./U))./((1+((6*fHz*340.2)./U)).^(5/3));
S = Sf/(2*pi);
SWind = S';

fname = [];

if Responsetype == 4
RAO = RAOP;
RAOpha = PhasePitch;
I = 0:0.01:30;
end
if Responsetype == 1
RAO = RAOTMy;                  
RAOpha = PhaseTMy;
I = 1000:1000:10000000;
end
if Responsetype == 2
RAO = RAOLC1;               
RAOpha = PhaseLoadcell1;
I = 10000:1000:100000000;
end
if Responsetype == 3
RAO = RAONx;                  
RAOpha = PhaseNx;
 I = 0:0.01:15;
end

wvec = f; % set angular frequencies
Period = (2*pi)./wvec;
SHz = SWind*(2*pi);

% Resample so that spectrum and RAO frequencies match up

%resample RAO at frequencies of S
fRAO = 1./fRAO; % as fRAO in Hz
fRAO(isinf(fRAO)|isnan(fRAO)) = fRAO(2)+(0.01*fRAO(2));
tsin = timeseries(RAO,fRAO);
tsout = resample(tsin,Period);
RAOsampledWind = flip(squeeze(tsout.data));

tsin = timeseries(RAOpha,fRAO);
tsout = resample(tsin,Period);
phaseofinterestWind = flip(squeeze(tsout.data));

SrWind = (SWind').*(RAOsampledWind'.^2); % convert to response spectrum. 
fWind = f;

% plot RAOs and resampled values to check  
% figure
% plot(fRAO,RAO)
% hold on
% plot(Period,(RAOsampledWind))
 
% figure
% plot(Period,SrWave)
% hold on
% yyaxis right
% plot(Period,(Sr))

% Calculate spectral moments and amplitudes
% S in frequency components (for M0)
      
b = (SrWind).*fc;


% for M1
% f = wvec   
b1 = (SrWind).*fc.*(f);


%M2       
b2= (SrWind).*fc.*(f.^2);

% M4   
b4 = (SrWind).*fc.*(f.^4);
     
LFbj = sqrt(2*SrWind.*fc);

% moments 
M0Wind = sum(b,2);
M1Wind = sum(b1,2);
M2Wind = sum(b2,2);
M4Wind = sum(b4,2);
MM0Wind = sqrt(M0Wind);


%% 2. calculate target amplitude
% from moments of response spectrum
ET = (((hourexposureLFMLER*60)*60)/(sqrt(scalefactor))); % Exposure time
ee = sqrt(1-(((M2Wind)^2)/(M0Wind*M4Wind))); % bandwidth parameter
n = (1/(4*pi))*((1+sqrt(1-(ee^2)))/sqrt(1-(ee^2)))*sqrt(M2Wind/M0Wind)*(ET); % (2.12) expected number of peaks in exposure time
TEV = sqrt(2*log((2*n*sqrt(1-(ee^2)))/(1+sqrt(1-(ee^2))))); % Target extreme value, a value of say 4 means a value at 4 times the standard deviation of the process 

nLF_old = n;

% Predict wind response EVD
FFmax =[];
eeee = [];
sig = 1;
muLF = TEV*sig;
m = round(1/(normcdf(0,muLF,sig))); % number of samples of gaussian process for exposure time

% rayleigh peaks 
% for i = 1:length(I)
%     eeee(i) = I(i)/MM0Wind;
% FFmax(i) = n*((I(i)/((MM0Wind^2)))*exp(-(I(i)^2)/(2*(MM0Wind^2)))*(1-exp(-(I(i)^2)/(2*(MM0Wind^2))))^(n-1));
% end
% Gaussian process
% for i = 1:length(I)
%     eeee(i) = I(i)/MM0Wind;
% FFmax(i) = m*((1/((MM0Wind*sqrt(2*pi))))*exp(-(I(i)^2)/(2*(MM0Wind^2))))*(normcdf(eeee(i),0))^(m-1);
% end
        eeee = I/MM0Wind;
        FFmax = m*((1/((MM0Wind*sqrt(2*pi)))).*exp(-(I.^2)/(2*(MM0Wind^2)))).*(normcdf(eeee,0)).^(m-1);

tarLF = cumtrapz(eeee*MM0Wind,FFmax); % get CDF of extreme value distribution
if plots == 1
plot(eeee*MM0Wind,tarLF,'r-.')
grid on
end
tarpksLF = nthroot(tarLF,n); % find CDF of peaks rather than extremes
LFXaxis = eeee*MM0Wind;

[minValue,closestIndex] = min(abs(tarLF-percentileComp099'));
MLF099 = eeee(closestIndex)*MM0Wind; % Target amplitude of response
[minValue,closestIndex] = min(abs(tarLF-percentileComp05'));
MLF05 = eeee(closestIndex)*MM0Wind; % Target amplitude of response
[minValue,closestIndex] = min(abs(tarLF-percentile'));
MLF = eeee(closestIndex)*MM0Wind; % Target amplitude of response

%% 4. Estimate the relative importance of the waves relative to wind and calculate a weighting
Relimp099 = (MWF099)/(MLF099); % relative importance = ratio of wind to wave response calculated from 99th percentiles
Relimp05 = (MWF05)/(MLF05); % relative importance = ratio of wind to wave response calculated from 50th percentiles

Absimp = MWF099 + MLF099; 
% make minus if wind dominates, positive if wave dominate, makes for a
% clearer plot.
if Relimp099 <1
Relimp099 = (1./Relimp099)*-1;
end
if Relimp05 <1
Relimp05 = (1./Relimp05)*-1;
end
Relimp099T = Relimp099;
% If wave dominates use 50th percentile for rel imp
if Relimp099 >0 && Relimp099 <RIC
    Relimp099 = Relimp05;
end

% calculate weighting to split combined target response
1/(abs(Relimp099)+1);
weighting = (ans*abs(Relimp099));
if Relimp099 <0
    weighting = 1-weighting;
end


%% 5. Estimate the target response extreme value distribution of the combined response from waves and wind  

% estimate of combined spectrum shape (location wrong) 
 CDF = ((weighting*tarpksWF)+((1-weighting)*tarpksLF)); % pks

% New way to get combined dist,
% calc number of samples (m) for exposure time
FFmax =[];
eeee = [];
sig = 1;
% The value of m needs to reflect the importance of wind relative to wave. 
% The wind is lower frequency and so there will be an over prediction of the response if the number of samples is calculated using the number of wave peaks or weighting them too heavily  
% This bit could be improved...
if Relimp099<-1 % if wind dominates 
 m = round(1/(normcdf(0,(muWF*weighting)+(muLF*(1-weighting)),sig)));
else
 m11 = round(1/(normcdf(0,muLF,sig)));
m22 = round(1/(normcdf(0,muWF,sig)));
m = m11+m22;  % if wave dominates
end

% estimate std of combined process
n = 50000000;     % Number of samples
% Generate Gaussian processes
process1 = (muWF) + MM0Wave  * randn(n, 1); % Wave
process2 = (muLF) + MM0Wind * randn(n, 1); % Wind
% Calculate the sum of the processes
sum_process = process1 + process2;
% calculate std of process
 MM0 = (std(sum_process));
 if Responsetype > 2
MM0 = round(std(sum_process),2); % observed to lead to a stable result from n = 50000000 simulations for this response
 else
  MM0 = round(std(sum_process),4,"significant"); % observed to lead to a stable result from n = 50000000 simulations for this response
 end   
 % Is there an analytical way to calculate MM0? (std of combined process)

 % combined EVD
% for i = 1:length(I)
%     eeee(i) = I(i)/MM0;
% FFmax(i) = m*((1/((MM0*sqrt(2*pi))))*exp(-(I(i)^2)/(2*(MM0^2))))*(normcdf(eeee(i),0))^(m-1);
% end
        eeee = I/MM0;
        FFmax = m*((1/((MM0*sqrt(2*pi)))).*exp(-(I.^2)/(2*(MM0^2)))).*(normcdf(eeee,0)).^(m-1);

tarGcomb = cumtrapz(eeee*MM0,FFmax); %check CDF goes to 1
if Relimp099 > 0
    if plots ==1
        grid on
plot(eeee*MM0,tarGcomb,'y-.') 
    end
end

% Shift EVD by 99th percentile of CDF as this is observed to have a better
% shape closer to simulations. I don't like this... This could also be improved (removing this step) by finding a
% better analytical way to get the combined EVD. Doesn't have any impact on
% the 99th percentile target but does improve 80th percentile target
% slightly. Method would still work if this step were removed.
if Relimp099 < 0
% find 99th percentile to shift dist....
[minValue,closestIndex] = min(abs(tarGcomb-0.99'));
valT = eeee(closestIndex)*MM0;
[minValue,closestIndex] = min(abs((CDF.^(nWF_old+nLF_old))-0.99'));
valT2 = LFXaxis(closestIndex)
Shift = valT-valT2;
if plots ==1
 plot(LFXaxis+Shift,CDF.^(nWF_old+nLF_old),'y-.')
end
end
if plots ==1
 plot(LFXaxis,CDF.^(nWF_old+nLF_old),'--k')
 legend('Wave','Wind','True combined','Combined 1st estimate')
end
 %% calculate the combined response target amplitude Total
[val,ind] = ((min(abs((tarGcomb-percentile)))));
if Relimp099 < 0
    [val,ind] = ((min(abs(((CDF.^(nWF_old+nLF_old))-percentile)))));
Total = LFXaxis(ind)+Shift;
else
    Total = LFXaxis(ind);
end
if plots == 1
xlim([Total*0.6 Total*1.2])
end

%split the target combined response into wind and wave contributions
%according to the relatie importance.
     if Relimp099>0 % if wave more important
Wavetar = (1/(1+Relimp099))*Relimp099*Total;
Windtar = (1-(1/(1+Relimp099))*Relimp099)*Total;
else
Windtar = (1/(1+abs(Relimp099)))*abs(Relimp099)*Total;
Wavetar = (1-(1/(1+abs(Relimp099)))*abs(Relimp099))*Total;
end
% end

if plots ==1
hold on
if Responsetype == 1
xlabel('TwrMy (Nm)')
xlim([100000 800000])
end
if Responsetype == 2
xlabel('Mooring load (N)')
xlim([100000 8000000])
end
if Responsetype == 3
xlabel('Nxa (m/s^2)')
xlim([0.5 4])
end
if Responsetype == 4
xlabel('Pitch (deg)')
xlim([1 8])
end

ylabel('CDF')
end

% Create MLER for Wind
 if abs(Relimp099)<RIC % split target percentiles between wave and wind.
M = Windtar;
 end
 if Relimp099<-RIC  % Just create wind profile
     M = MLF;
 end
 if abs(Relimp099)<RIC || Relimp099<-RIC 

     % calculate wind profile
fr = M1Wind/M0Wind;
fn = fWind; % frequency components

% Spectral amplitude of wind Spectrum
% for i = 1:length(SWind)
% aa(i) = sqrt((SWind(i)*fc(i)));   % according to NREL  application of MLER paper
% end
aa = sqrt((SWind'.*fc));   % according to NREL  application of MLER paper

phi = sqrt(SrWind'./SWind);
phi(isnan(phi))=0;

% for i = 1:length(SrWind)
% A(i) = ((phi(i)*aa(i))/((M0Wind*M2Wind)-(M1Wind^2)))*((M2Wind-(M1Wind*fn(i)))+(fr*((M0Wind*fn(i))-M1Wind)));
% end
A = ((phi'.*aa)./((M0Wind*M2Wind)-(M1Wind^2))).*((M2Wind-(M1Wind.*fn))+(fr.*((M0Wind.*fn)-M1Wind)));

% wind profile
phe = phaseofinterestWind; 


 jj = abs(min(ttMLER))*2; % MLER peaks at half this value
ZZZ = zeros(1, length(ddt:ddt:jj));
jjj =  ddt:ddt:jj;
 for j = 1:length(jjj)
Z = (aa.*M.*A.*(cos((fn.*(jjj(j)-(jj/2)))+(phe')))); 
Z = Z(~isnan(Z));
 ZZZ(j) = sum(Z);
 end

 j = ddt:ddt:jj;
[v ind] = min(abs(j-(abs(min(ttMLER))+max(ttMLER)))); %j - total length of run
% plot((j(1:ind)-jj./2),ZZZ(1:ind))
j = j(1:ind);
ZZZ = ZZZ(1:ind);
% figure
% hold on
% plot((j(1:ind)-jj./2),ZZZ(1:ind)+U,'r','linewidth',1.5)
% grid on
% xlabel('Time (s)')
% ylabel('Wind speed (m/s)')

%%
% % save in matlab format
P = [];
P = cat(1,(j-j(1)),(ZZZ+U));
P = P';

% uniform wind input file used in openfast
a = zeros(1,length(ZZZ));
Shear = linspace(0, 0,length(a)); % set the shear wanted for OpenFAST
% Shear = linspace(0.11, 0.11,length(a)); % power law exponent
W = cat(1,P(:,1)',P(:,2)',a,a,a,Shear,a,a);
W = W';

if abs(Relimp099)<RIC
 fname = strcat(sprintf(fnamMLERWind),'Split.mat'); % matlab (Time Series)
 fname1 = strcat(sprintf(fnamMLERWind),'Split.csv'); % csv (Time Series)
   fname2 = strcat(sprintf(fnamMLERWind),'Split.wnd'); % OpenFAST (Time Series)
else 
 fname = strcat(sprintf(fnamMLERWind),'.mat'); % matlab (Time Series)
 fname1 = strcat(sprintf(fnamMLERWind),'.csv'); % csv (Time Series)
   fname2 = strcat(sprintf(fnamMLERWind),'.wnd'); % OpenFAST (Time Series)
end

  if saveq == 1
  % save as matlab
   save(fname,'P','Relimp099','Relimp099T', 'Relimp05','Absimp','Wavetar','Windtar','Total','XX') 
 % OpenFAST format
% dlmwrite(fname2,W,'delimiter',' ','precision',7); % time series
fid1 = fopen(fname2, 'w');
            fprintf(fid1, '%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\n', W');
            fclose(fid1);
  else
  end
PMLERWind = P;
%% %%%%%%%%%%%%%%%%%%%%%%% Create CMLERs (CRRWs) for Wind

phe = phaseofinterestWind;
phi = RAOsampledWind;

%surface elevation
aaa = sqrt((SWind'.*fc));   

%response
aa = sqrt((SrWind.*fc));  

% var = (aa.^2); 
sig = 1;  % Dietz pg 49 (pg 46 for single MLER)
mu = 0;
        Vn = normrnd(mu, sig, length(aa), CRRWnum);
        Wn = normrnd(mu, sig, length(aa), CRRWnum);

% 
S1 = sum( aa'.*((Vn.*cos(phe))+(Wn.*sin(phe))), 1 );
S2 = sum(aa'.*fn'.*((Vn.*sin(phe))-(Wn.*cos(phe))));
S3 = sum(aa'.*((-Vn.*sin(phe))+(Wn.*cos(phe))));
S4 = sum(aa'.*fn'.*((Vn.*cos(phe))+(Wn.*sin(phe))),1);

%% Condition V and W
fm = M1Wind/M0Wind;
             V = Vn - ((aa'/((M2Wave*M0Wave)-M1Wave^2)).*(((M2Wave-(fn'*M1Wave)).*S1.*cos(phe))-(M * (M2Wave-(fn'*M1Wave)).*cos(phe))+(((fn'*M0Wave)-M1Wave).*S2.*sin(phe))+(((fn'*M1Wave)-M2Wave).*S3.*sin(phe))+(((fn'*M0Wave)-M1Wave).*S4.*cos(phe))-(((M * fm*((fn'*M0Wave)-M1Wave).*cos(phe))))));
        W = Wn - ((aa'/((M2Wave*M0Wave)-M1Wave^2)).*(((M2Wave-(fn'*M1Wave)).*S1.*sin(phe))-(M * (M2Wave-(fn'*M1Wave)).*sin(phe))-(((fn'*M0Wave)-M1Wave).*S2.*cos(phe))-(((fn'*M1Wave)-M2Wave).*S3.*cos(phe))+(((fn'*M0Wave)-M1Wave).*S4.*sin(phe))-(((M * fm*((fn'*M0Wave)-M1Wave).*sin(phe))))));


SE = [];
SEA = zeros(CRRWnum, length(ttMLER));  
        for t = 1:length(ttMLER)
            SE = aaa'.*((V.*cos((-fn'.*ttMLER(t))))+(W.*sin((-fn'.*ttMLER(t)))));
            SEA(:,t) = sum(SE,1);
        end

ttMLER2 = 0:ddt:(ttMLER(end)-(ttMLER(1)));  
repmat(ttMLER2,size(SEA,1),1);
SSE = cat(3,ans,SEA);

if plots ==1
figure
for i = 1:CRRWnum
    hold on
plot(SSE(i,:,1),(SSE(i,:,2)+U))
end
plot(SSE(1,:,1),mean(SSE(:,:,2))+U,'k','linewidth',1.5)
end
CRRWMean2 = cat(1,squeeze(SSE(1,:,1)),mean(SSE(1:CRRWnum,:,2),1)+U);
CRRWMean2=CRRWMean2';

a = zeros(1,size(CRRWMean2,1));
Shear = linspace(0, 0,length(a)); 
% Shear = linspace(0.11, 0.11,length(a)); % power law exponent
CRRWW = cat(1,CRRWMean2(:,1),CRRWMean2(:,2),a',a',a',Shear',a',a');
% CRRWW = CRRWW';

if plots == 1
grid on 
xlabel('Time (s)')
ylabel('Wind speed (m/s)')
% plot(PMLERWind(:,1),PMLERWind(:,2),'r','linewidth',1.5)
plot(CRRWMean2(:,1),CRRWMean2(:,2),'r','linewidth',1.5)

xlim([250 500])
end
%% Save CRRW for paddle software and matlab

for i = 1:CRRWnumSave
fname = [];
P = cat(1,squeeze(SSE(i,:,1))/sqrt(scalefactor),(squeeze(SSE(i,:,2))+U)/scalefactor);
P = P';
a = zeros(1,size(P,1));
Shear = linspace(0, 0,length(a)); 
% Shear = linspace(0.11, 0.11,length(a)); % power law exponent
W = cat(1,P(:,1)',P(:,2)',a,a,a,Shear,a,a);
W = W';

if abs(Relimp099)<RIC
 fname = strcat(fnamCMLERWind,sprintf('Split_%d',i),'.mat'); %matlab (Time Series)
%  fname1 = strcat(fnamCRRWwind,sprintf('Split_%d',i),'.csv'); % csv (Time Series)
   fname2 = strcat(fnamCMLERWind,sprintf('Split_%d',i),'.wnd'); % OpenFAST (Time Series)
%   fname1 = strcat(sprintf('CRRW_LC2_099_%d',i),'.csv'); % 
fname3 = strcat(fnamCMLERWindmean,sprintf('Split',i),'.mat'); % instead of MLER for wind, use mean of CMLER cases
fname4 = strcat(fnamCMLERWindmean,sprintf('Split',i),'.wnd'); % instead of MLER for wind, use mean of CMLER cases
else
   fname = strcat(fnamCMLERWind,sprintf('_%d',i),'.mat'); %matlab (Time Series)
%  fname1 = strcat(fnamCRRWwind,sprintf('_%d',i),'.csv'); % csv (Time Series)
   fname2 = strcat(fnamCMLERWind,sprintf('_%d',i),'.wnd'); % OpenFAST (Time Series)
%   fname1 = strcat(sprintf('CRRW_LC2_099_%d',i),'.csv'); 
fname3 = strcat(fnamCMLERWindmean,'.mat'); % instead of MLER for wind, use mean of CMLER cases
fname4 = strcat(fnamCMLERWindmean,'.wnd'); % instead of MLER for wind, use mean of CMLER cases
end

if saveq == 1
save(fname,'P','Relimp099T','Relimp099', 'Relimp05','Absimp','Wavetar','Windtar','Total') 
save(fname3,'CRRWMean2','Relimp099T','Relimp099', 'Relimp05','Absimp','Wavetar','Windtar','Total') 

% save OpenFAST wind files
% dlmwrite(fname1,P,'delimiter','\t','precision',7);
 dlmwrite(fname2,W,'delimiter',' ','precision',7); % time series
 fid1 = fopen(fname2, 'w');
            fprintf(fid1, '%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\n', W');
            fclose(fid1);
 a = zeros(1,size(CRRWMean2,1));
Shear = linspace(0, 0,length(a)); 
% Shear = linspace(0.11, 0.11,length(a)); % power law exponent
CRRWW = cat(1,CRRWMean2(:,1)',CRRWMean2(:,2)',a,a,a,Shear,a,a);
CRRWW = CRRWW';
  % dlmwrite(fname4,CRRWW,'delimiter',' ','precision',7); % time series
                fid1 = fopen(fname4, 'w');
              fprintf(fid1, '%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\t%g\t%10.8f\n', CRRWW');
            fclose(fid1);
else
end
P = []; 
end

 end

%% Produce wave cases

% if Relimp >-1*RIC || splitchoice == 1
%%%%%%%%%%%%%%%%%%%%%%%%% Create MLER for Waves
 if abs(Relimp099)<RIC % split target percentiles between wave and wind.
M = Wavetar;
 if (Relimp099)>0 % split target percentiles between wave and wind.
M = Wavetar;
 end
 end
 if Relimp099>RIC
     M = MWF;
 end

%% 3. calculate wave profile
fr = M1Wave/M0Wave;
fn = fWave;
fc = diff(fn);
% Spectral amplitude of sea Spectrum
aa = sqrt((SWave*fc(1)));   % according to NREL  application of MLER paper

phi = sqrt(SrWave./SWave);
phi(isnan(phi))=0;

A = ((phi.*aa)./((M0Wave*M2Wave)-(M1Wave^2))).*((M2Wave-(M1Wave.*fn'))+(fr*((M0Wave.*fn')-M1Wave)));

%% wave profile, at position x

% wave number
 g = 9.81; 
 fHZ = fn/(2*pi);
 for i = 1:length(fHZ)
    L(i) = (g*((1/fHZ(i))^2)/(2*pi));
    La(i) = (g*((1/fHZ(i))^2)/(2*pi))*tanh((2*pi*d)/L(i));
        while abs(La(i) - L(i))>0.01
            L(i) = La(i);
            La(i) = (g*((1/fHZ(i))^2)/(2*pi))*tanh((2*pi*d)/La(i)); %wavelength
        end
        % Calculate the wave number, group and phase velocities
        K(i) = ((2*pi)/La(i));
%         c = sqrt((g/k)*tanh(k*d));
%         cg = (c/2)*(1+((2*k*d)/sinh(2*k*d)));
 end
K(isnan(K))=K(2);


phe = phaseofinterestWave; 
ZZZ = [];
Z = [];
jj = abs(min(ttMLER))*2; % MLER peaks at half this value
ZZZ = zeros(1, length(ddt:ddt:jj));
jjj =  ddt:ddt:jj;
 for j = 1:length(jjj)
Z = (aa.*M.*A.*(cos((fn'.*(jjj(j)-(jj/2)))+(phe)))); 
Z = Z(~isnan(Z));
 ZZZ(j) = sum(Z);
 end
 j = ddt:ddt:jj;
[v ind] = min(abs(j-(abs(min(ttMLER))+max(ttMLER)))); % j - total length of run
if plots ==1
    figure
hold on
plot((j(1:ind)-jj./2),ZZZ(1:ind))
hold on
grid on
xlabel('Time (s)')
ylabel('Surface elevation  (m)')
end
j = j(1:ind);
ZZZ = ZZZ(1:ind);


% save in matlab format
P = [];
P = cat(1,(j-j(1)),ZZZ);
P = P';
P = P(1:end-1,:);  %needs to be an odd number if importing as spectrum to
PMLERWave = P;
% wecsim.
if plots ==1
plot(P(:,1)-1,P(:,2))
end
 if abs(Relimp099)<RIC % split target percentiles between wave and wind.
 fname = strcat(sprintf(fnamMLER),'Split.mat'); % matlab (Time Series)
 fname1 = strcat(sprintf(fnamMLER),'Split.csv'); % wavepaddles (Time Series)
   fname2 = strcat(sprintf(fnamMLER),'Split.Elev'); % OpenFAST (Time Series)
 else
 fname = strcat(sprintf(fnamMLER),'.mat'); % matlab (Time Series)
 fname1 = strcat(sprintf(fnamMLER),'.csv'); % wavepaddles (Time Series)
   fname2 = strcat(sprintf(fnamMLER),'.Elev'); % OpenFAST (Time Series)
 end

  % save 
  if saveq == 1
   save(fname,'P','Relimp099T','Relimp099', 'Relimp05','Absimp','Wavetar','Windtar','Total','XX') 
% dlmwrite(fname2,P,'delimiter','\t','precision',7); % time series
fid1 = fopen(fname2, 'w');
            fprintf(fid1, '%g\t%10.8f\n',P);
            fclose(fid1);
  else
  end
%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create CMLERs (CRRWs) for Waves

phe = phaseofinterestWave;
phi = RAOsampledWave;

SAA = [];
SAAF = [];
SAAB= [];

%surface elevation
aaa = sqrt((SWave*fc(1)));   

%response
aa = sqrt((SrWave*fc(1)));  %not sure why the 2 is left out? same for MLER paper. 

% var = (aa.^2); 
sig = 1;  % Dietz pg 49 (pg 46 for single MLER)
mu = 0;
        Vn = normrnd(mu, sig, length(aa), CRRWnumSave);
        Wn = normrnd(mu, sig, length(aa), CRRWnumSave);

% 
        S1 = sum( aa.*((Vn.*cos(phe))+(Wn.*sin(phe))), 1 );
S2 = sum(aa.*fn'.*((Vn.*sin(phe))-(Wn.*cos(phe))));
S3 = sum(aa.*((-Vn.*sin(phe))+(Wn.*cos(phe))));
S4 = sum(aa.*fn'.*((Vn.*cos(phe))+(Wn.*sin(phe))),1);

%% Condition V and W
fm = M1Wave/M0Wave;

             V = Vn - ((aa/((M2Wave*M0Wave)-M1Wave^2)).*(((M2Wave-(fn'*M1Wave)).*S1.*cos(phe))-(M * (M2Wave-(fn'*M1Wave)).*cos(phe))+(((fn'*M0Wave)-M1Wave).*S2.*sin(phe))+(((fn'*M1Wave)-M2Wave).*S3.*sin(phe))+(((fn'*M0Wave)-M1Wave).*S4.*cos(phe))-(((M * fm*((fn'*M0Wave)-M1Wave).*cos(phe))))));
        W = Wn - ((aa/((M2Wave*M0Wave)-M1Wave^2)).*(((M2Wave-(fn'*M1Wave)).*S1.*sin(phe))-(M * (M2Wave-(fn'*M1Wave)).*sin(phe))-(((fn'*M0Wave)-M1Wave).*S2.*cos(phe))-(((fn'*M1Wave)-M2Wave).*S3.*cos(phe))+(((fn'*M0Wave)-M1Wave).*S4.*sin(phe))-(((M * fm*((fn'*M0Wave)-M1Wave).*sin(phe))))));

SE = [];
SEA = zeros(CRRWnumSave, length(ttMLER));  
        for t = 1:length(ttMLER)
            SE = aaa.*((V.*cos((-fn'.*ttMLER(t))))+(W.*sin((-fn'.*ttMLER(t)))));
            SEA(:,t) = sum(SE,1);
        end
ttMLER = 0:ddt:(ttMLER(end)-(ttMLER(1)));  
repmat(ttMLER,size(SEA,1),1);
SSE = cat(3,ans,SEA);

if plots ==1
figure
for i = 1:CRRWnumSave
    hold on
plot(SSE(i,:,1),SSE(i,:,2))
end
plot(SSE(1,:,1),mean(SSE(:,:,2)),'k','linewidth',1.5)
grid on 
plot(PMLERWave(:,1),PMLERWave(:,2),'r','linewidth',1.5)
xlim([300 450])
end

% save
for i = 1:CRRWnumSave
fname = [];
   P = SSE(i,:,:);
   P = squeeze(P);
    if abs(Relimp099)<RIC % split target percentiles between wave and wind.
 fname = strcat(sprintf(fnamCMLER),sprintf('Split_%d',i),'.mat'); %matlab (Time Series)
%  fname1 = strcat(sprintf(fnam),sprintf('_%d',i),'.csv'); % wavepaddles (Time Series)
   fname2 = strcat(sprintf(fnamCMLER),sprintf('Split_%d',i),'.Elev'); % OpenFAST (Time Series)
    else
 fname = strcat(sprintf(fnamCMLER),sprintf('_%d',i),'.mat'); %matlab (Time Series)
%  fname1 = strcat(sprintf(fnam),sprintf('_%d',i),'.csv'); %wavepaddles (Time Series)
   fname2 = strcat(sprintf(fnamCMLER),sprintf('_%d',i),'.Elev'); % OpenFAST (Time Series)
    end
   if saveq == 1
   save(fname,'P','Relimp099T','Relimp099', 'Relimp05','Absimp','Wavetar','Windtar','Total','XX') ;
% save in paddle software format
% dlmwrite(fname2,P,'delimiter','\t','precision',7);
fid1 = fopen(fname2, 'w');
            fprintf(fid1, '%g\t%10.8f\n', P);
            fclose(fid1);
  else
  end
P = []; 

end

if plots ==1 && saveq == 1
 figure(1)
 legend('Wave','Wind','Combnew','Comb check')
 % save figure of target EVDs 
    H = figure(1);
savefig(H,strcat(fnamEVD,'.fig'));
end

Extime = toc
