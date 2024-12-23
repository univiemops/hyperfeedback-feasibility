%% Preprocessing for music data 
% 
% preprocessing using homer 2 functions, partly adapted ("_adapt" postfix)
% contains: channel pruning, conversion to optical density, TDDR motion
% correction with low-pass filtering, conversion to concentration, baseline
% PCA. 
% follows the real-time recommendations of Klein (2024).
%
% KK last version October 2024

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
%
% hbo was used in paper
%if one chromophore is bad, we will not be able to distinguish
%between hbo and hbr. Therefore we need to block the whole channel.


Npair=12;

pathn='/path/to/data/';
savepathn='/path/to/preprocessed_data/';

% count sh*t channels
shitchanlist = zeros(1,3);
counter = 1;


for npair=1:Npair
 for cond=1%:2 %only using part learning
     for nprobe=1:4
         if cond ==1
             %probe
             filen=['FF_' num2str(npair) '_MES_Probe' num2str(nprobe)];
         elseif cond ==2
             filen=['FF' num2str(cond) '_' num2str(npair) '_MES_Probe' num2str(nprobe)];
         end
            
    load([pathn filen '.nirs'], '-mat');
     
   disp(['Preprocessing data of ' filen '...']);


% % % % % prune channels

% checking for bad channels and removing them (SD.MeasListAct has zeros 
% input for bad channels)

 tInc      = ones(size(d,1),1);                                                 
 dRange    = [0.01 4];
 SNRthresh = 2;
 resetFlag = 0;
 SD.MeasListAct =  ones(size(SD.MeasList,1),1); 
 SDrange = [20 40];

 
 SD       = enPruneChannels(d, SD, tInc, dRange,...
                                 SNRthresh, SDrange, resetFlag);
                             
% block both chromophores even though maybe only one is classified as bad:
removechan = find(SD.MeasListAct==0);
SD.MeasListAct(removechan + 22) = 0;
removechan2 = [removechan; removechan+22];


% % % % % convert to optical density

dod = hmrIntensity2OD(d);   


% % % % % TDDR motion correction with low-pass filter

fs             = 1/(t(2)-t(1));  % sampling frequency of the data
dodTDDR = hmrMotionCorrectTDDR_adapt(dod, SD, fs);


% convert changes in OD to changes in concentrations (HbO, HbR, and HbT)
ppf      = [6 6];   % partial pathlength factors for each wavelength.
dc       = hmrOD2Conc(dodTDDR, SD, ppf);


% % % % % baseline PCA filter

%PCA filter from homer2
%[fdc, ~, ~] = hmrPCAFilter(dc, SD, [1,1]);

%PCA baseline filter
baseline = dc(1:600,:,:); %take first minute as baseline for baseline PCA
dc_task = dc(601:end, :, :);
t = t(601:end);

%baseline PCA
[bdc, ~,~ ] = hmrPCAFilter_baseline(dc_task, baseline, SD, [1,1] );


% % % % % prepare and write out

% extract hbo and hbr
hbo=squeeze(bdc(:,1,:));
hbr=squeeze(bdc(:,2,:));

% NaN faulty channels
 hbo(:,removechan')=NaN(size(hbo,1),size(removechan,1));
 hbr(:,removechan')=NaN(size(hbo,1),size(removechan,1));

for m=1:length(removechan)%length(SD.MeasListAct) %for counting, see below

        sh1tchanlist(counter,1) = npair;
        sh1tchanlist(counter,2) = nprobe;
        sh1tchanlist(counter,3) = removechan(m); 
        counter=counter + 1;

end
  
       save([savepathn filen '.mat'], 'hbo','hbr','s','t', 'fs');
 
     end
 end
end



%% info about faulty/ pruned channels

%sort sh*t channels to paper numbering
sh1tchanlist_red = sh1tchanlist;

for s=1:size(sh1tchanlist,1)
    if sh1tchanlist(s,2) == 1 || sh1tchanlist(s,2) == 3
        sh1tchanlist_red(s,2) = 1; %subject 1
    elseif sh1tchanlist(s,2) == 2 || sh1tchanlist(s,2) == 4
        sh1tchanlist_red(s,2) = 2; %subject 2
    else
        warning('list mismatch in faulty channel list, please check');
    end
    if sh1tchanlist(s,2) == 2 || sh1tchanlist(s,2) == 4 %don't distinguish between subjects bc if one sub is bad whole dyad is bad
        sh1tchanlist_red(s,3) =  sh1tchanlist(s,3) +22;
    else
        %keep as is
    end
end

%to be excluded channels for main analysis: every channel which is bad for
%more than three subjects/ dyads
excludeChanList = zeros(2,1);
counter = 1;

faultChan = unique(sh1tchanlist_red(:,3));
c=categorical(sh1tchanlist_red(:,3), faultChan);
[histN, histChannels] = histcounts(c);

%cutoff 3 dayds
cut = find(histN>=3);
rlist = histChannels(cut);



