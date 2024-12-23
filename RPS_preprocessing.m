%% Preprocessing for RPS data 
% 
% preprocessing using homer 2 functions, partly adapted ("_adapt" postfix)
% contains: channel pruning, conversion to optical density, TDDR motion
% correction with low-pass filtering, conversion to concentration, baseline
% PCA. 
% follows the real-time recommendations of Klein (2024).
%
% KK last version October 2024


%% FP

clear all
srcPath = '/path/to/rawdata';                        % raw data location
desPath = '/path/to/processed/data';                  % processed data location

% for baseline PCA, get markers from files for third rest condition, use
%rest_indices3 = restMat(3,1):restMat(3,2)-1;

% count sh*t channels
sh1tchanlist_1={};
sh1tchanlist_2={};


pairlist   = {'01','02','03','04','05','06','08','09','10','11','13', ...
    '14','15','17','18','20','21', '22', '23','24','26','28','29', ...
    '30','31','32'};
nPairs  = length(pairlist);

%conditions
conds = {'FP','PS', 'PD','C'};
co    = 1; %using FP for greatest contrast
%sampling rate
fs      = 7.8125; %Hz

%% preprocessing sub1


for pair=1:length(pairlist)
    
    filen = ['/RPS_' pairlist{pair} '_sub1_' conds{co}];
    
    load([srcPath filen '.nirs'], '-mat');
    
    load([srcPath '/RPS_' pairlist{pair} '_marker_' ...
        conds{co} '.mat']); %get marker matrices
    
    %indices of last rest block
    rest_indices3 = restMat(3,1):restMat(3,2)-1;
    %remove 8 sec for equilibrium
    rest_indices= restMat(3,1)+round(8*fs):restMat(3,2)-1;
    
    
% % % % % prune channels    
                         
% checking for bad channels and removing them (SD.MeasListAct has zeros 
% input for bad channels)

 tInc      = ones(size(d,1),1);                                                 
 dRange    = [0.01 2.5];
 SNRthresh = 2;
 resetFlag = 0;
 SD.MeasListAct =  ones(size(SD.MeasList,1),1); 
 SDrange = [2.0 4.5];

 
 SD       = enPruneChannels(d, SD, tInc, dRange,...
                                 SNRthresh, SDrange, resetFlag);

                             
% convert the wavelength data to optical density
dod = hmrIntensity2OD(d);                               
                             


% % % % % TDDR motion correction with low-pass filter

dodTDDR = hmrMotionCorrectTDDR_adapt(dod, SD, fs);
                                    

% convert changes in OD to changes in concentrations (HbO, HbR, and HbT)
ppf      = [6 6];                                                       % partial pathlength factors for each wavelength.
dc       = hmrOD2Conc(dodTDDR, SD, ppf);


% % % % % baseline PCA filter

% %PCA filter from homer2
% [fdc, svs, nSV] = hmrPCAFilter(dc, SD, [1,1]);


%split data set for baseline PCA
baseline = dc(rest_indices, :, :);

dc_task = dc(1:restMat(3,1),:,:);
t       = t(1:restMat(3,1));
%baseline PCA
[bdc, ~,~ ] = hmrPCAFilter_baseline(dc_task, baseline, SD, [1,1] );


% % % % % prepare and write out

% extract hbo and hbr
hbo=squeeze(bdc(:,1,:));
hbr=squeeze(bdc(:,2,:));


for m=1:length(SD.MeasListAct)/2 % other half is other chromophore but we NaN anyways
    if SD.MeasListAct(m) == 0
        sh1tchanlist_1 =[sh1tchanlist_1; {filen, m}]; 
        %NaN pruned channels
        hbo(:,m) = NaN;
        hbr(:,m) = NaN;
    else
        %keep channels
    end
end
   
% save preprocessed data
save([desPath filen '.mat'], 'hbo','hbr','s','t', 'fs');


end



%% preprocessing sub 2

for pair=1:length(pairlist)
    
    filen = ['/RPS_' pairlist{pair} '_sub2_' conds{co}];
    
    load([srcPath filen '.nirs'], '-mat');
    
    load([srcPath '/RPS_' pairlist{pair} '_marker_' ...
        conds{co} '.mat']); %get marker matrices
    
    %indices of last rest block
    rest_indices3 = restMat(3,1):restMat(3,2)-1;
    %remove 8 sec for equilibrium
    rest_indices= restMat(3,1)+round(8*fs):restMat(3,2)-1;
                          
    
% % % % % prune channels    

% checking for bad channels and removing them (SD.MeasListAct has zeros 
% input for bad channels)

 tInc      = ones(size(d,1),1);                                                 
 dRange    = [0.01 2.5];
 SNRthresh = 2;
 resetFlag = 0;
 SD.MeasListAct =  ones(size(SD.MeasList,1),1); 
 SDrange = [2.0 4.5];

 
 SD       = enPruneChannels(d, SD, tInc, dRange,...
                                 SNRthresh, SDrange, resetFlag);
                             
% convert the wavelength data to optical density
dod = hmrIntensity2OD(d);                          
                           

% % % % % TDDR motion correction and low-pass filter
% 
dodTDDR = hmrMotionCorrectTDDR_KKadapt(dod, SD, fs);

% convert changes in OD to changes in concentrations (HbO, HbR, and HbT)
ppf      = [6 6];                                                       % partial pathlength factors for each wavelength.
dc       = hmrOD2Conc(dodTDDR, SD, ppf);


% % % % % baseline PCA filter

%PCA filter from homer2
%[fdc, ~, ~] = hmrPCAFilter(dc, SD, [1,1]);


%split data set for baseline PCA
baseline = dc(rest_indices, :, :);

dc_task = dc(1:restMat(3,1),:,:);
t       = t(1:restMat(3,1));
%baseline PCA
[bdc, ~,~ ] = hmrPCAFilter_KKbaseline(dc_task, baseline, SD, [1,1] );



% % % % % prepare and write out

%extract hbo and hbr
hbo=squeeze(bdc(:,1,:));
hbr=squeeze(bdc(:,2,:));

% 
for m=1:length(SD.MeasListAct)/2 % other half is other chromophore but we NaN anyways
    if SD.MeasListAct(m) == 0
        sh1tchanlist_2 =[sh1tchanlist_2; {filen, m}]; 
        %NaN pruned channels
        hbo(:,m) = NaN;
        hbr(:,m) = NaN;
    else
        %keep channels
    end
end
   
% save preprocessed data

save([desPath filen '.mat'], 'hbo','hbr','s','t', 'fs');
end


