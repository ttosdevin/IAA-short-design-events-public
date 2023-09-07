% This script loads in matlab variables from a wave only white noise or
% wind only (with parked turbine) simulation and processes the RAOs as
% needed to run the short design events script. The matlab variables were
% extracted from OpenFAST outputs.
clear
clc
close all

% load data
ijk = 1 % 1 is Windcrete
for kji = 1:2 % 1 = wave, 2 = wind
if kji ==2
     fnam = 'RAO_U50' % wind
else
   fnam = 'WhiteNoise' % wave
end
load(fnam)

%% Filtering
data = FT2;
  cut = [35]/(2*pi); %cut off frequency
  bb = 1/(0.5*diff(Time(1:2))); %half the sampling frequency
  wn = cut/bb;
  [b,a] = butter(1,wn,'low');
  FT2 = filtfilt(b,a,data);

if exist('FT1') 
    data = FT1; % FT1. FT4 for 5mw barge. 
  cut = [35]/(2*pi); %cut off frequency
  bb = 1/(0.5*diff(Time(1:2))); %half the sampling frequency
  wn = cut/bb;
  [b,a] = butter(1,wn,'low');
  FT1 = filtfilt(b,a,data);
end
if exist('FT3') 
      data = FT3; % FT1. FT4 for 5mw barge. 
  cut = [35]/(2*pi); %cut off frequency
  bb = 1/(0.5*diff(Time(1:2))); %half the sampling frequency
  wn = cut/bb;
  [b,a] = butter(1,wn,'low');
  FT3 = filtfilt(b,a,data);
end

%% get RAO magnitudes
ff = 0:0.001:20; % angfreq range 
 WS = 10000; % windowing
 cut = 400; % time to cut from start of run before analysis 
datrange = [round(cut/Time(2)):length(Surge)]; % Time series to use
b = (Time(end)/length(Time));
 Fs=1./(b); % sample freq

% Window size in [] should be selected carefully...
% output freq in Hz
if kji == 1
[EncounteredWaveSpectrum,fe] = pwelch(Wave(datrange),[WS],[],[ff/(2*pi)], Fs); % wave
else
[EncounteredWaveSpectrum,fe] = pwelch(Windx(datrange),[WS],[],[ff/(2*pi)], Fs); % wind
end

[ResponseSpectrumP,fP] = pwelch(Pitch(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumNx,fNx] = pwelch(Nxa(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumNz,fNz] = pwelch(Nza(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumS,fS] = pwelch(Surge(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumLC1,fLC1] = pwelch(FT2(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumTMy,fTMy] = pwelch(TwrMy(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumH,fH] = pwelch(Heave(datrange),[WS],[],[ff/(2*pi)], Fs);
[ResponseSpectrumTFx,fFx] = pwelch(TwrFx(datrange),[WS],[],[ff/(2*pi)], Fs);

RAOP = sqrt(ResponseSpectrumP./EncounteredWaveSpectrum);
RAOS = sqrt(ResponseSpectrumS./EncounteredWaveSpectrum);
RAOH = sqrt(ResponseSpectrumH./EncounteredWaveSpectrum);
RAONx = sqrt(ResponseSpectrumNx./EncounteredWaveSpectrum);
RAONz = sqrt(ResponseSpectrumNz./EncounteredWaveSpectrum);
RAOLC1 = sqrt(ResponseSpectrumLC1./EncounteredWaveSpectrum);
RAOTMy = sqrt(ResponseSpectrumTMy./EncounteredWaveSpectrum);
RAOTFx = sqrt(ResponseSpectrumTFx./EncounteredWaveSpectrum);

IND=find(EncounteredWaveSpectrum > (.03*max(EncounteredWaveSpectrum))); % RAOs unrealistic outside wave/wind frequency region 
INDstart = IND(1);
INDend = IND(end);


% Plot RAOs
% figure
% plot(1./(fLC1(IND(1):IND(end))),RAOLC1(IND(1):IND(end)));
%  hold on
% plot((1./fP(IND(1):IND(end))),RAOP(IND(1):IND(end)));

% Plot spectra
% figure
% plot(fe,EncounteredWaveSpectrum*2)
% hold on
% plot(fP,ResponseSpectrumP*2);

%% calculate the phase of the RAOs using cpsd
if kji == 1
[Pxy,fRAO] = cpsd(Wave(datrange),Pitch(datrange),[WS],[],[ff/(2*pi)], Fs);
PhasePitch   = angle(Pxy);
[Pxy,fRAO] = cpsd(Wave(datrange),Heave(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseHeave   = angle(Pxy);
[Pxy,fRAO] = cpsd(Wave(datrange),Surge(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseSurge   = angle(Pxy);
[Pxy,fRAO] = cpsd(Wave(datrange),Nxa(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseNx   = angle(Pxy); 
[Pxy,fRAO] = cpsd(Wave(datrange),Nza(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseNz   = angle(Pxy);
[Pxy,fRAO] = cpsd(Wave(datrange),FT2(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseLoadcell1  = angle(Pxy);
%   plot(1./(fff(IND(1):IND(end))),(PhaseLoadcell1(IND(1):IND(end))))
[Pxy,fRAO] = cpsd(Wave(datrange),TwrMy(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseTMy  = angle(Pxy);
[Pxy,fRAO] = cpsd(Wave(datrange),TwrFx(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseTFx  = angle(Pxy);
else
   [Pxy,fRAO] = cpsd(Windx(datrange),Pitch(datrange),[WS],[],[ff/(2*pi)], Fs);
PhasePitch   = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),Heave(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseHeave   = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),Surge(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseSurge   = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),Nxa(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseNx   = angle(Pxy); %245712-245704
[Pxy,fRAO] = cpsd(Windx(datrange),Nza(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseNz   = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),FT2(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseLoadcell1  = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),TwrMy(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseTMy  = angle(Pxy);
[Pxy,fRAO] = cpsd(Windx(datrange),TwrFx(datrange),[WS],[],[ff/(2*pi)], Fs);
PhaseTFx  = angle(Pxy);
end 

%% save
     fnamesave = strcat(sprintf(fnam),'_10000_RAOs');
      save(fnamesave)
end