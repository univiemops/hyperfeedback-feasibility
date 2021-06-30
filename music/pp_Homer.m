%% PP for music data adapted from pp for RPS data
%for now with Homer2, potentially to be adapted to Homer 3
%
% KK 18.06.21 need to decide what to do with bad intensities (0 or <0).
% Filter function does not accept NaN. Probably need to throw out channels,
% via adapting measlistact file (and adapt filter so that this file gets
% considered)
%
%
%KK 15.06.21 basics are working. For now need to do: 
% 1. write script adapting Hitachi_2_Homer_v3 to get all participants into
% nirs format.
% 2. preprocess all subjects according to Trinh's script.
% 3. look at wtc of desired channels (check significant channels in paper)


% %% FP
% 
% clear all
% srcPath = 'P:\projects\RPS\RPS\procData\nirsData\nirsData_FP\';                        % raw data location
% desPath = 'P:\projects\RPS\RPS\procData\hmrData\Data_FP\';                  % processed data location
% 
% 
% %% Scan for all recordings
%   sourceList    = dir([srcPath, 'RPS_*_sub1_FP.nirs']);
%   sourceList    = struct2cell(sourceList);
%   sourceList    = sourceList(1,:);
%   numOfSources  = length(sourceList);
%   PartList       = zeros(1, numOfSources);
%   PartList  = erase(sourceList, '.nirs');
% 
% 
% %% preprocessing
% for i = 1:1:numOfSources
% 
%   % load raw data of subject 1
%   srcFolder   = strcat(srcPath);
%   filename    = PartList{1, i};
%   
%   fprintf('Load raw nirs data of subject...\n');
%     file_path = strcat(srcFolder, filename,'.nirs');

Npair=12;

pathn='/Users/kathrin/data/Music_fNIRS/data/';

%faulty testfile
%load([pathn 'FF2_2_MES_Probe4.nirs'], '-mat');

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


%problem: some intensity values are negative or 0. Trying to 'NaN' them...
problemvec=find(d<=0);
dNaN=d;
dNaN(problemvec)=NaN;

%NaNs are not accepted for the filter. For now block the channel...
[~,col]=find(d<=0); %find all columns(=channels) with bullsh*t values
removecol=unique(col);
dZero=d;
dZero(:,removecol)=0; %zero the channels andhope the filter doesn't crash... nope, does crash


dod = hmrIntensity2OD(dZero);    %writes NaN in dod                       


% checking for bad channels and removing them (SD.MeasListAct has zeros 
% input for bad channels)

 tInc      = ones(size(d,1),1);                                                 
 dRange    = [0 10000000];
 SNRthresh = 2;
 resetFlag = 0;
 SD.MeasListAct =  ones(size(SD.MeasList,1),1); 
 SDrange = [2.5 3.5];

 
 SD       = enPruneChannels(dod, SD, tInc, dRange,...
                                 SNRthresh, SDrange, resetFlag);

                         
                             
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
dodSpline(:,removecol)=0;


% % Identify bad channels manually
% ppf   = [6 6];                                                              % partial pathlength factors for each wavelength.
% dod_check  = hmrOD2Conc(dodSpline, SD, ppf);
% hbo = squeeze(dod_check(:,1,:));
% for i = 1:1:size(hbo, 2)
%   subplot(6,4,i);
%   if ~isnan(hbo(:,i))
%   sig = [t, hbo(:,i)];
%   sigma2=var(sig(:,2));                                                     % estimate signal variance
%   
%   [wave,period,~,coi,~] = wt(sig);                                          % compute wavelet power spectrum
%   power = (abs(wave)).^2 ;
%   
%   for j=1:1:length(coi)
%     wave(period >= coi(j), j) = NaN;                                        % set values below cone of interest to NAN
%   end
% 
%   h = imagesc(t, log2(period), log2(abs(power/sigma2)));
%   colorbar;
%   Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
%   set(gca,'YLim',log2([min(period),max(period)]), ...
%           'YDir','reverse', 'layer','top', ...
%           'YTick',log2(Yticks(:)), ...
%           'YTickLabel',num2str(Yticks'), ...
%           'layer','top')
%   title(sprintf('Channel %d', i));
%   ylabel('Period in seconds');
%   xlabel('Time in seconds');
%   set(h, 'AlphaData', ~isnan(wave));
% 
%   colormap jet;
%   end
%   end
% set(gcf,'units','normalized','outerposition',[0 0 1 1])                     % maximize figure
% 
% badChannels = channelCheckbox();
                             
    

% correcting for motion artifacts using Wavelet-based motion correction.                                

% iQr             = 1.5;
% [~, dod_corr]  = evalc(...                                             % evalc supresses annoying fprintf output of hmrMotionCorrectWavelet
%                 'hmrMotionCorrectWavelet(dod, SD, iQr);');






% bandpass filtering
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


% 
%    
% % save preprocessed data
%   desFolder   = strcat(desPath);
% %   filename    = sprintf(PartList{1, i});
%   
%   file_path = strcat(desFolder, filename, ...
%                      '.mat');
% 
%   fprintf('The preprocessed data of dyad will be saved in'); 
%   fprintf('%s ...\n', file_path);

   save([pathn filen '.mat'], 'hbo','hbr','s','t', 'fs');

%   fprintf('Data stored!\n\n');
%   clear data_raw
% end


% 
     end
 end
end