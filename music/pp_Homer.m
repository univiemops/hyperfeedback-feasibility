%% Preprocessing for music data adapted from preprocessing for RPS data
% with Homer2
%
% KK June 2021, bugfix Jan 2023

%FF is the prefix for the “part learning” condition
%FF2 is the prefix for the “whole learning” condition
%Probe1 = right hemisphere sub 1
%Probe2 = left hem sub 1
%Probe3 = right hem sub 2
%Probe4 = left hem sub2
%time points 1:1800 for resting state and 1801:6600 for task phase
%from paper: 
%FOI: 0.07-0.15Hz  (period 6.61-14.01s)
%musical phrase 6.64\pm 1.56s
%Nch=44
% compute complete rest block, take blocks of task phase
% taking hbo, as in paper




Npair=12;

pathn='YourPathToData/nirsData/'; %Unix system used, beware syntax
savepathn='YourPathToData/ppData/';



for npair=1:Npair
 for cond=1:2
     for nprobe=1:4
         if cond ==1
             filen=['FF_' num2str(npair) '_MES_Probe' num2str(nprobe)];
         elseif cond ==2
             filen=['FF' num2str(cond) '_' num2str(npair) '_MES_Probe' num2str(nprobe)];
         end
    

         
    load([pathn filen '.nirs'], '-mat');
    
   
   disp(['Preprocessing data of ' filen '...']);
 
   
   
% convert the wavelength data to optical density

%if one chromophore is bad, we will not be able to distinguish
%between hbo and hbr. Therefore we need to block the whole channel.



dod = hmrIntensity2OD(d);    %writes NaN in dod                       


% checking for bad channels and removing them (SD.MeasListAct has zeros 
% input for bad channels)

 tInc      = ones(size(d,1),1);                                                 
 dRange    = [0 10000000];
 SNRthresh=0.01;
 resetFlag = 0;
 SD.MeasListAct =  ones(size(SD.MeasList,1),1); 
 SDrange = [20 40];

 
 SD       = enPruneChannels(dod, SD, tInc, dRange,...
                                 SNRthresh, SDrange, resetFlag);
                             
% %problem: some intensity values are negative or 0. We need to eliminate
% these bad channels.
[~,col]=find(d<=0); %find all columns with bullsh*t values
removecol=unique(col);
disp(removecol);

dod(:,removecol)=0; %for smooth further processing use '0', NaN in the end.

                            
                             
% identifies motion artifacts in an input data matrix d. If any active
% data channel exhibits a signal change greater than std_thresh or
% amp_thresh, then a segment of data around that time point is marked as a
% motion artifact.
tInc            = ones(size(d,1),1);                                                 
tMotion         = 1;
tMask           = 1;
stdevThreshold  = 50;
ampThreshold    = 0.4;
fs             = 1/(t(2)-t(1));                              % sampling frequency of the data
[tIncAuto,tIncCh]       =  hmrMotionArtifactByChannel(dod, fs, SD,...
                                        tInc, tMotion,...
                                        tMask, stdevThreshold,...
                                        ampThreshold);

                                    
% % % % Spline interpolation
p=0.99;
dodSpline = hmrMotionCorrectSpline(dod, t, SD, tIncCh, p);   


%here comes the shitty workaround to get rid of the NaNs that make the
%filter crash
%dodSpline(:,removecol)=0;



% % % % bandpass filtering
lpf             = 0.5;                                                  % in Hz
hpf             = 0.01;                                                 % in Hz


dod_corr_filt  = hmrBandpassFilt(dodSpline, fs, hpf, ...
                                      lpf);

% convert changes in OD to changes in concentrations (HbO, HbR, and HbT)
ppf      = [6 6];                                                       % partial pathlength factors for each wavelength.
dc       = hmrOD2Conc(dod_corr_filt, SD, ppf);

% extract hbo and hbr
hbo = squeeze(dc(:,1,:));
hbr = squeeze(dc(:,2,:));


%NaN bad channels to be on the safe side

removechan=removecol;
for r=1:length(removecol)
    if removecol(r)>22 %22 channels, 44 bc of two chromophores. This is only valid for this specific data set!
        removechan(r)=removecol(r)-22;
    end
end
removechan=unique(removechan);

hbo(:,removechan')=NaN(size(hbo,1),size(removechan,1));
hbr(:,removechan')=NaN(size(hbo,1),size(removechan,1));

clear removecol removechan;

% 


       save([savepathn filen '.mat'], 'hbo','hbr','s','t', 'fs');
% 



% 
     end
 end
end
