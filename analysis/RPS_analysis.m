%% coherence calculation and statistical test for the RPS data set.
% This script calculates the wavelet tranform coherence within a certain
% frequency band (period of interest) for different task windows. It then
% calculates the wavelet tranform coherence within a certain period band of
% interest (frequency band of interest) for different task windows. It then
% calculates the difference in coherence between the task window and the
% complete rest condition (excluding the data used for the baseline PCA
% during preprocessing). Statistics are assessed via a general linear mixed
% model under the null hypothesis that there is no difference between task
% data and rest data. All data are corrected for multiple comparisons using
% an FDR routine. Single channel as well as summary plots are generated.
%
% KK last version October 2024

path    = '/path/to/preprocesseddata';
    
%markerpath:
mpath   = '/path/to/markers/'; 

pairlist   = {'01','02','03','04','05','06','08','09','10','11','13', ...
    '14','15','17','18','20','21', '22','23','24','26','28','29', ...
    '30','31','32'}; 

nPairs  = length(pairlist);

%sampling rate
fs      = 7.8125; %Hz

%conditions
conds = {'FP','PS', 'PD','C'};
co    = 1; %using FP for greatest contrast

windowsizelist = [round(50*fs) round(60*fs) round(70*fs) round(80*fs) ...
    round(90*fs) round(100*fs)]; %in samples
nTestWindows = length(windowsizelist);

%offset between windows
woffset = round(8*fs); %in samples, offset between windows

%total number of channels: 16, channel 1-8 DLPFC, channel 9-16 TPJ
channellist = 1:16;
nCh = length(channellist);

% period band of interest
poi = [6 14];

% chromophore switch
chromophore = 1; %1= HbO, 2=HbR; HbR is used in Kayhan et al. (2022)

% switch for writing out data and figures
writedata = true;
saveworkspace = true;
savefig = false;

%initialize for stats:
%7 outputs to glme.coefficients
teststats_glme     = NaN(nTestWindows, nCh, 7); 
all_teststats_glme = NaN(nTestWindows+2, nCh, 7);


%initialize for plots:
res_fig_mean = NaN(nTestWindows+2, nCh);
res_fig_std  = NaN(nTestWindows+2, nCh,2);



%get exact parameters for period and filter

%get data of exemplary pair
if chromophore==1
    
    load([path '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test = hbo;
    load([path  '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test = hbo;
    
elseif chromophore==2
    load([path  '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test = hbr;
    load([path '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test = hbr;
    
else
    error('chromophore not specified.');
end

%get exemplary markers - marker matrices
load([mpath '/RPS_01_marker_' conds{co} '.mat']); 
if ~exist('t') || ~exist('fs') ||  ~exist('restMat') || ...
        ~exist('taskMat') || ~exist('subtaskMat')
    error(['failed to load in markers and/or data for pair1']);
end

% get indices form experiment markers
task_indices1 = subtaskMat(1,1):restMat(2,1)-1;
task_indices2 = subtaskMat(1,2):restMat(3,1)-1;
rest_indices1 = restMat(1,1):restMat(1,2)-1;
rest_indices2 = restMat(2,1):restMat(2,2)-1;
rest_indices3 = restMat(3,1):restMat(3,2)-1;


sigPart1 = [t, hb1test(:,14)];
sigPart2 = [t, hb2test(:,14)];


%get cutoff for coi for Matlab wavelet toolbox:
[~,~, periodtry,coitry] = wcoherence(sigPart1(:,2), sigPart2(:,2), ...
    seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], ...
    'VoicesPerOctave', 14);

% calculate coi 
 f=1./seconds(max(periodtry));
 cf = 6/(2*pi);
 predtimes = sqrt(2)*cf./f; % coi per wavelet frequency; in seconds
 coiidx= round(predtimes*fs); %in samples


for win = 1:length(windowsizelist)
    disp(win);
    
    windowsize = windowsizelist(win);
    
    %number of windows
    nWinRest1 = floor(length(rest_indices1) / (windowsize + woffset));
    nWinRest2 = floor(length(rest_indices2) / (windowsize + woffset));
    nWinRest3 = floor(length(rest_indices3) / (windowsize + woffset));
    
    nWinTask1 = floor(length(task_indices1) / (windowsize + woffset));
    nWinTask2 = floor(length(task_indices2) / (windowsize + woffset));
    
    nWinRest = nWinRest1 + nWinRest2 + nWinRest3;
    nWinTask = nWinTask1 + nWinTask2;
    
    %resulting coherence values
    c_rest_win = NaN(nWinRest, nPairs, nCh);
    c_task_win = NaN(nWinTask, nPairs, nCh);
    
    c_rest_tot = NaN(3, nPairs, nCh);
    c_task_tot = NaN(2, nPairs, nCh);
 
    c_rest_tot_offmat = NaN(3, nPairs, nCh);
    c_task_tot_offmat = NaN(2, nPairs, nCh);
    
    %matrices for the differences between task and rest coherence
    diff_task_win = NaN(nWinTask, nPairs, nCh);
    diff_task_tot = NaN(2, nPairs, nCh);
    diff_task_off = NaN(2, nPairs, nCh);
    


    %load in pairs, use part learning condition
    for pair=1:length(pairlist)
        
        %load in participants
        
        if chromophore==1
            
            load([path '/RPS_' pairlist{pair} ...
                '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1 = hbo;
            load([path  '/RPS_' pairlist{pair} ...
                '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2 = hbo;
            
        elseif chromophore==2
            
            load([path  '/RPS_' pairlist{pair} ...
                '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1 = hbr;
            load([path  '/RPS_' pairlist{pair} ...
                '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2 = hbr;
            
        else
            error('chromophore not specified.');
        end
        
        %in case markers differ slightly for this pair, extract them again:
        load([mpath '/RPS_' pairlist{pair} '_marker_' ...
            conds{co} '.mat']); %get marker matrices
        if ~exist('t') || ~exist('fs') ||  ~exist('restMat') || ...
                ~exist('taskMat') || ~exist('subtaskMat')
            error(['failed to load in markers and/or data for pair ' ...
                num2str(pairlist{pair})]);
        end
        task_indices1 = subtaskMat(1,1):restMat(2,1)-1;
        task_indices2 = subtaskMat(1,2):restMat(3,1)-1;
        rest_indices1 = restMat(1,1):restMat(1,2)-1;
        rest_indices2 = restMat(2,1):restMat(2,2)-1;
        rest_indices3 = restMat(3,1):restMat(3,2)-1;
        
        
        %~~~~ now calculate coherence for every channel
        for channel=1:length(channellist)%Nch
            
            ch=channellist(channel);
            
            % get signal of interest = soi
            % check if this channel was not rejected in both subjects 
            % during preprocessing
            if ~isnan(hb1(1, ch)) && ~isnan(hb2(1, ch))  
                sigPart1 = [t, hb1(:,ch)];
                sigPart2 = [t, hb2(:,ch)];
                
            else
                warning(['Houston, we are having a problem in channel ' ...
                    num2str(ch) 'for dyad ' pairlist{pair}]);
                continue
            end
            
            
            %~~~~ online windowed coherence in windows

            w_counter = 0;
            
            for w = 1:nWinTask1
                w_counter = w_counter + 1;
                startpoint = task_indices1(1) + woffset; %in samples
                part = (startpoint + (w-1)*woffset + (w-1)*windowsize): ...
                    (startpoint + (w-1)*woffset + w*windowsize - 1);
                %calc coherence
                c_task_win(w_counter, pair, channel) = ...
                    calculateCoherence(sigPart1, sigPart2, fs, poi, ...
                    coiidx, part);
            end
            
            for w = 1:nWinTask2
                w_counter = w_counter + 1;
                startpoint = task_indices2(1) + woffset; %in samples
                part = (startpoint+(w-1)*woffset+(w-1)*windowsize): ...
                    (startpoint+(w-1)*woffset+w*windowsize-1);
                c_task_win(w_counter,pair,channel) = ...
                    calculateCoherence(sigPart1, sigPart2, fs, poi,...
                    coiidx, part);
            end
            
            %~~~~ online complete blocks
            
            %rest blocks
            part = rest_indices1;
            c_rest_tot(1,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            part = rest_indices2;
            c_rest_tot(2,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            
            %task blocks
            part =task_indices1;
            c_task_tot(1,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            part =task_indices2;
            c_task_tot(2,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            
            %~~~~~ offline complete timecourse
            
            %calc coherence
            [pRsq, ~, period, coi] = wcoherence(sigPart1(:,2), ...
                sigPart2(:,2), seconds(1/fs),'PeriodLimits', ...
                [seconds(6) seconds(14)], 'VoicesPerOctave', 14);
            
            for j = 1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN; 
            end
            
            %reduced pRsq to avoid mayer wave freq
            redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
            
            % Get signal of interest = soi
            soi = squeeze(mean(redpRsq, 1, 'omitnan'));
            
 
            
            soi_task1 = soi(task_indices1);
            soi_task2 = soi(task_indices2);
            soi_rest1 = soi(rest_indices1);
            soi_rest2 = soi(rest_indices2); 
            

            
            c_rest_tot_offmat(1, pair, channel) = mean( ...
                soi_rest1(coiidx:length(soi_rest1)-coiidx), 'omitnan');
            c_rest_tot_offmat(2, pair, channel) = mean( ...
                soi_rest2(coiidx:length(soi_rest2)-coiidx), 'omitnan');

            c_task_tot_offmat(1, pair, channel) = mean( ...
                soi_task1(coiidx:length(soi_task1)-coiidx), 'omitnan');
            c_task_tot_offmat(2, pair, channel) = mean( ...
                soi_task2(coiidx:length(soi_task2)-coiidx), 'omitnan');
            
            

% create differences 
            counter = 0;
            
            for j = 1:nWinTask1
                
                counter = counter + 1;
                diff_task_win(counter, pair, channel) = ...
                    c_task_win(counter,pair,channel) - ...
                    c_rest_tot(1,pair,channel);                               
            end
            
            for j = 1:nWinTask2
                
                counter = counter + 1;
                diff_task_win(counter, pair, channel) = ...
                    c_task_win(counter,pair,channel) - ...
                    c_rest_tot(2,pair,channel);
            end
            
            
            for j = 1:2
                
                diff_task_tot(j,pair,channel) = ...
                    c_task_tot(j,pair,channel)-c_rest_tot(j,pair,channel);
                
                diff_task_off(j,pair,channel) = ...
                    c_task_tot_offmat(j,pair,channel) - ...
                    c_rest_tot_offmat(j,pair,channel);
            end

            
        end %end channel
        
    end %end pairs
 
    % hierarchical binomial model; see function code below.
    teststats_glme(win,:,:) = calculateGlmeForChannels( ...
        diff_task_win, nCh, nPairs, nWinTask);
    all_teststats_glme(win,:,:) = teststats_glme(win,:,:);

end %end window loop


% %total data, online fashion

teststats_glme_tot = calculateGlmeForChannels( ...
    diff_task_tot, nCh, nPairs, 2); % nWin is two here because two runs.


% %total data, offline fashion

teststats_glme_off = calculateGlmeForChannels( ...
    diff_task_off, nCh, nPairs, 2); % nWin is two here because two runs.


all_teststats_glme(nTestWindows+1,:,:) = teststats_glme_tot;
all_teststats_glme(nTestWindows+2,:,:) = teststats_glme_off;


%% FDR routine

%q<0.05
testvec    = reshape(all_teststats_glme(:,:,5), (nTestWindows+2)*nCh, 1);
[n_signif05, index_signif05] = fdr(testvec, 0.05);
resultsvec = NaN(size(testvec));
resultsvec(index_signif05) = testvec(index_signif05);
resultsmat = reshape(resultsvec, (nTestWindows+2), nCh);
all_teststats_glme_FDR = all_teststats_glme;
all_teststats_glme_FDR(:,:,8) = resultsmat;

%q<0.01
[n_signif01, index_signif01] = fdr(testvec, 0.01);
resultsvec = NaN(size(testvec));
resultsvec(index_signif01) = testvec(index_signif01);
resultsmat = reshape(resultsvec, (nTestWindows+2), nCh);
all_teststats_glme_FDR(:,:,9) = resultsmat;

%% save workspace
if saveworkspace
    if chromophore == 1
    save('workspaceRPS_conFB_glme_FDR_revpp_HbO');
    elseif chromophore == 2
    save('workspaceRPS_conFB_glme_FDR_revpp_HbR');        
    else
    warning('Problem with chromophore identity.');
    end
else
    %do nothing
end



%% plot with error bars
% make total grandmean histogram - all channels and conjunction window
% significant and offline significant

    for i=1:size(all_teststats_glme,1)
        if i<=nTestWindows
            windows(i)={[num2str(round(windowsizelist(i)/fs)) ' sec']};
        elseif i==nTestWindows+1
            windows(i)={'complete'};
        elseif i==nTestWindows+2
            windows(i)={'offline'};
        end
    end

%number significant offline channels:
NSignOff=length(find(~isnan(all_teststats_glme_FDR(length(windows),:,8))));
chSignOff = find(~isnan(all_teststats_glme_FDR(length(windows),:,8)));

%count how many significant if offline significant
NSignBoth = zeros(length(windows),1);
for win = 1:length(windows)
    NSignBoth(win) = length(find(~isnan(all_teststats_glme_FDR(length(windows),:,8)) & ~isnan(all_teststats_glme_FDR(win,:,8))));
end


for win=1:nTestWindows+2
    for channel=1:nCh   
     
    % heights of bars:    
    res_fig_mean(win,channel) = (all_teststats_glme(win,channel,1));

    res_fig_std(win,channel,1) = all_teststats_glme(win,channel,2);
    res_fig_std(win,channel,2) = all_teststats_glme(win,channel,2);     
    end
    
end

    avgWTCsignif_win = NaN(length(windowsizelist)+2,1);
    WTCsignif_std_win = NaN(length(windowsizelist)+2,1);
    
    for swin=1:length(windowsizelist)+2
        
    idxsignif = find(~isnan(all_teststats_glme_FDR(length(windows),:,8)) & ~isnan(all_teststats_glme_FDR(win,:,8)));
    avgWTCsignif_win(swin) = mean(all_teststats_glme_FDR(swin, idxsignif, 1));
    WTCsignif_std_win(swin) = mean(all_teststats_glme_FDR(swin, idxsignif, 2));    
    
    end


% make a figure of average delta_WTC across all channels

% axis labels

    X=categorical(windows);
    X=reordercats(X,windows); %apparently Matlab is reordering the vector when 
    %transforming to categorical... so we need to change back.
    
    res_fig_grandmean = mean(res_fig_mean,2);
    res_fig_grandstd  = mean(res_fig_std(:,:,1),2);
    
    bar_data = [res_fig_grandmean  avgWTCsignif_win];
    error_data = [res_fig_grandstd  WTCsignif_std_win];
    
        figure;
    b1=bar(bar_data, 'grouped','LineWidth', 2, 'FaceAlpha', 0.8); hold on
    
    [ngroups,nbars] = size(bar_data);
    x = nan(nbars, ngroups);
    
for i = 1:nbars
    x(i,:) = b1(i).XEndPoints;
   % 
end
    er =errorbar(x',bar_data,error_data,'k','linestyle','none', 'LineWidth', 2);
    

    hold off
    
    yline(0,'LineWidth', 2);     
    
    axes=gca;
    axes.FontSize=18;
    axes.FontWeight='b';
    axes.LineWidth=2.5;
    ylim([0 0.15]);
    
    
    if chromophore==1
        title(['HbO, grandmean over all channels '], ...
            'FontSize', 18, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/Grandmeantot_HbO_' num2str(ch)], 'svg');
        else
            %do nothing
        end        
    elseif chromophore==2       
        title(['HbR, grandmean over all channels '], ...
            'FontSize', 18, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/Grandmeantot_HbR_' num2str(ch)], 'svg');
        else
            %do nothing
        end        
    else
        warning('Problem with chromophore identity.');
    end
    
%% make figure mean over offline significant channels

    res_fig_avgSig = mean(res_fig_mean(:,chSignOff),2);
    res_fig_stdSig  = mean(res_fig_std(:,chSignOff,1),2);%std(res_fig_mean(:,chSignOff),1,2);
      
        figure;
    bar(X,res_fig_avgSig(:),'LineWidth', 3); hold on
    er = errorbar(X,res_fig_avgSig,res_fig_stdSig, ...
        res_fig_stdSig,'LineWidth', 3);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    yline(0,'LineWidth', 4); 
    
    axes=gca;
    axes.FontSize=16;
    axes.FontWeight='b';
    axes.LineWidth=3;
    ylim([0 0.2]);
    
       for i=1:size(all_teststats_glme,1)       
    text(X(i), er.YData(i)+ er.YPositiveDelta(i)+0.05, [num2str(NSignBoth(i)) '\' num2str(NSignOff)], 'FontSize',15, 'FontWeight', 'b', 'HorizontalAlignment','center');
   % text(X(i), er.YData(i)+0.1, 'muh', 'FontSize',15, 'FontWeight', 'b', 'HorizontalAlignment','center');
    end
    
    hold off
    if chromophore==1
        %title({'difference in WTC, errorbars SEM' , ['for HbO, channel ' ...
        title(['HbO, grandmean over all significant channels '], ...
            'FontSize', 16, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/GrandmeanSign_HbO_' num2str(ch)], 'svg');
        else
            %do nothing
        end        
    elseif chromophore==2       
        title(['HbR, grandmean over all significant channels '], ...
            'FontSize', 16, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/GrandmeanSign_HbR_' num2str(ch)], 'svg');
        else
            %do nothing
        end        
    else
        warning('Problem with chromophore identity.');
    end  

%%

%make tiledlayout which includes figure for every channel %make a figure for every channel

fig=figure;
%fig.Position(3:4)=[1500,2000];
tlt=tiledlayout(4,4,'TileSpacing', 'compact');

for channel=1:nCh
    ch=channellist(channel);
    
    for i=1:size(all_teststats_glme,1)

        if i<=nTestWindows
            windows(i)={[num2str(round(windowsizelist(i)/fs)) ' sec']};
        elseif i==nTestWindows+1
            windows(i)={'complete'};
        elseif i==nTestWindows+2
            windows(i)={'offline'};
        end
    end
    nexttile; %figure;
    %X=categorical(windows);
    %X=reordercats(X,windows); %apparently Matlab is reordering vector when 
    %transforming to categorical... so we need to change back.
    
    %to save plot remove long categorical labels for all tiles except
    %for last
%     if channel < nCh
%     X=categorical(1:length(windows));
%     elseif channel == nCh
%     X=categorical(windows);
    X=reordercats(X,windows); %apparently Matlab is reordering vector when 
    %transforming to categorical... so we need to change back.  
%     end
    
    b=bar(X,res_fig_mean(:,channel),'LineWidth', 2.5); hold on
    er = errorbar(X,res_fig_mean(:,channel),res_fig_std(:,channel,2), ...
        res_fig_std(:,channel,1),'LineWidth', 2.5);  
    b.FaceColor = [0.231/0.255, 0.168/0.255, 0.168/0.255];
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    yline(0,'LineWidth', 2.5); 
    
    for i=1:size(all_teststats_glme,1)   
        if ~isnan(all_teststats_glme_FDR(i,channel,9)) 
            text(X(i), er.YData(i)+0.05, '\ast\ast', 'FontSize',22, 'FontWeight', 'b', 'HorizontalAlignment','center');
        elseif ~isnan(all_teststats_glme_FDR(i,channel,8)) 
            text(X(i), er.YData(i)+0.05, '\ast', 'FontSize',22, 'FontWeight', 'b', 'HorizontalAlignment','center');
        else
            %no star
        end
        
    end
    
    axes=gca;
    axes.FontSize=15;%25;
    axes.FontWeight='b';
    axes.LineWidth=3;
  %  yticks([0 0.1 0.15]);
    ylim([0 0.15]);
    
    hold off

        title(['channel ' ...
            num2str(ch)], 'FontSize', 15, 'FontWeight', 'b');    
    
end


if chromophore==1
    title(tlt,'\Delta WTC per channel and window size for the Rock-Paper-Scissors data, HbO', 'FontSize', 20, 'FontWeight', 'b');
    if savefig
        saveas(gcf, ['./figures/contFB_allchan_HbO_' num2str(ch)], 'svg');
    else
        %do nothing
    end
elseif chromophore==2
    title(tlt,'\Delta WTC per channel and window size for the Rock-Paper-Scissors data, HbR', 'FontSize', 20, 'FontWeight', 'b');
    if savefig
        saveas(gcf, ['./figures/contFB_allchan_HbR_' num2str(ch)], 'svg');
    else
        %do nothing
    end
else
        warning('Problem with chromophore identity.');
end
%% write out data - glme mixed model

if writedata

for channel=1:nCh
    %     windows=num2cell(round(windowsizelist/fs));
    ch=channellist(channel);
    for i=1:nTestWindows+2
        if i==nTestWindows+1
            windows(i)={'complete block'};
        elseif i==nTestWindows+2
            windows(i)={'offline'};
        else
            windows(i)={['windowsize ' ...
                num2str(round(windowsizelist(i)/fs)) ' sec']};
        end
    end
    rows={'estimate'; 'SE'; 'tStat'; 'DF'; 'pValues'; 'Lower'; 'Upper'; 'q<0.05'; 'q<0.01'};
    mytable=array2table(squeeze(all_teststats_glme_FDR(:,channel,:))', ...
        'VariableNames', windows, 'RowNames', rows);
    
    if chromophore==1
    writetable(mytable,['./tables/music_glme_contFB_FDR' num2str(ch) 'HbO_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.csv'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/music_glme_contFB_FDR' num2str(ch) 'HbR_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.csv'],'WriteRowNames',true)    
    else
        warning('Problem with chromophore identity.');
    end
end

else 
    %do nothing
end




%% function definitions

% coherence
function c_value = calculateCoherence(sigPart1, sigPart2, fs, poi, ...
    coiidx, part)

    % Calculate coherence with Matlab's coherence function
    [pRsq, ~, period, coi] = wcoherence(sigPart1(part, 2), ...
        sigPart2(part, 2), seconds(1/fs), 'PeriodLimits', ...
        [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);

    % NaN everything outside of coi as safety net
    for j = 1:1:length(coi)
        pRsq(period >= coi(j), j) = NaN;
    end
    
    %reduced pRsq to avoid mayer wave freq
    redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
    
    % Get signal of interest = soi
    soi = squeeze(mean(redpRsq, 1, 'omitnan'));
  %  soi = soi(part);

    c_value = mean(soi(coiidx:length(soi) - coiidx), 'omitnan');
end


% mixed model
function teststats_glme = calculateGlmeForChannels(diff_task_win, nCh,...
    nPairs, nWin)

    teststats_glme = NaN(nCh, 7);

    for channel = 1:nCh

        y = squeeze(diff_task_win(:, :, channel));
        
        % Initialize the result matrix for this channel
        mat = NaN(nPairs * nWin, 2);

        for a = 1:nPairs
            for b = 1:nWin
                
                 mat((a - 1) * nWin + b , 1) = atanh(y(b, a));                

                mat((a - 1) *nWin+b, 2) = a;
            end
        end
           

        % Create a table from the result matrix
        tab = array2table(mat, 'VariableNames', {'DeltaWTC', 'dyadID'});

        % Fit a generalized linear mixed model
        glme = fitglme(tab, 'DeltaWTC ~ 1 + (1|dyadID)', 'Distribution', ...
            'normal');

        % Due to data set format, we have to convert first to table, then
        % array
        conv = dataset2table(glme.Coefficients);
        teststats_glme(channel, :) = table2array(conv(1, 2:8));
        
    end
end
