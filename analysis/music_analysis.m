%% coherence calculation and statistical test for the music data set.
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



path=['/path/to/preprocessed_data/']; 
%use your data path; Unix System used, adapt accordingly

%experiment specifics
pairlist={'01','02','03','04','05','06', '07','08','09','10','11','12'};
nPairs=length(pairlist);

totCh=44; %Number of channels; bad channels will have zeros
fs=10; %sampling rate 10 Hz
task_indices=1800:6600-600; %signal of interest for task block, first minute of demo not used
rest_indices=1:1800-600; %experiment starts with rest, first minute used for baseline PCA
poi=[6 14.0]; %poi = period of interest

%windowlengths to be examined in samples
windowsizelist=[500 600 700 800 900 1000]; %in samples
nTestWindows=length(windowsizelist);


woffset=round(8*fs); %in samples, 8 sec offset between windows


%if only significant channels wanted:
channellist=[1,2,3,4,5,6,7,8,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25,26,28,29,31,32,34,35,36,37,39,40,41,42,44];
nCh=length(channellist);

%initialize variables for results of sign test 
%(coherence_task > coherence_rest) and resulting signs (VZ)
results=NaN(4,nTestWindows,nCh); %(p,zval,sign, total);
results_tot=NaN(4,nCh);
results_off=NaN(4,nCh);


%teststats_glm=NaN(nTestWindows,nCh,4); %4 outputs to mdl.coefficients: estimate, SE, tstat, pValue
teststats_glme=NaN(nTestWindows,nCh,7); %7 outputs to glme.coefficients
all_teststats_glme = NaN(nTestWindows+2, nCh, 7);


%initialize for plots:
res_fig_mean = NaN(nTestWindows+2, nCh);
res_fig_std  = NaN(nTestWindows+2, nCh,2);


%chromophore switch
chromophore = 2; %1= HbO, 2=HbR
%in paper HbO is used

writedata     = true;
savefig       = false;
saveworkspace = true;
    
%get exact parameters for period and filter
    
    %get data of exemplary pair
    if chromophore==1
        load([path '/FF_3_MES_Probe1.mat']);
        hb1test=hbo;
        load([path '/FF_3_MES_Probe3.mat']);
        hb2test=hbo;
    elseif chromophore==2
        load([path '/FF_3_MES_Probe1.mat']);
        hb1test=hbr;
        load([path '/FF_3_MES_Probe3.mat']);
        hb2test=hbr;
    else
        error('chromophore not specified.');
    end
    
    
    pnoi = zeros(2,1);
    
    sigPart1 = [t, hb1test(:,14)];
    sigPart2 = [t, hb2test(:,14)];

%get cutoff for coi for Matlab wavelet toolbox:    
[~,~, periodtry,coitry]=wcoherence(sigPart1(:,2), sigPart2(:,2), ...
    seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], ...
    'VoicesPerOctave', 14);


% calculate coi 
 f=1./seconds(max(periodtry));
 cf = 6/(2*pi);
 predtimes = sqrt(2)*cf./f; % coi per wavelet frequency; in seconds
 coiidx= round(predtimes*fs); %in samples

%%
for win=1:length(windowsizelist)
    disp(win);
    
    windowsize=windowsizelist(win);
    
    %number of windows
    nWinRest=floor(length(rest_indices)/(windowsize+woffset)); %number 
    %of possible windows per dyad 
    nWinTask=floor(length(task_indices)/(windowsize+woffset));
    
   % nWin_task_vec(win)=nWinTask;
    
    %resulting coherence values
    c_rest_win=NaN(nWinRest, nPairs, nCh);
    c_task_win=NaN(nWinTask, nPairs, nCh);

    c_rest_tot=NaN(nPairs, nCh);
    c_task_tot=NaN(nPairs, nCh);
    
    c_rest_tot_offmat=NaN(nPairs, nCh);
    c_task_tot_offmat=NaN(nPairs, nCh);
    
    %matrices for the differences between task and rest coherence
    diff_task_win=NaN(nWinTask, nPairs, nCh);
    diff_task_tot=NaN(1,nPairs,nCh);
    diff_task_off=NaN(1,nPairs,nCh);
    
    
    
        %load in pairs, use part learning condition
    for pair=1:length(pairlist)
        
        %load in participants
        
        hb1=zeros(size(hb1test,1),totCh);
        hb2=zeros(size(hb2test,1),totCh);
        
        
        if chromophore==1
            
            %right hemisphere both participants
            %get data (hb and fs)
            load([path '/FF_' num2str(pair) '_MES_Probe1.mat']); 
            hb1(:,1:totCh/2)=hbo; %first 22 channels right hem
            load([path '/FF_' num2str(pair) '_MES_Probe3.mat']);
            hb2(:,1:totCh/2)=hbo;
            
            %left hemisphere both participants
            %get data (hb and fs)
            load([path '/FF_' num2str(pair) '_MES_Probe2.mat']); 
            hb1(:,totCh/2+1:totCh)=hbo; %next 22 channels left hem
            load([path '/FF_' num2str(pair) '_MES_Probe4.mat']); 
            hb2(:,totCh/2+1:totCh)=hbo;
            
        elseif chromophore==2
            
            %right hemisphere both participants
            %get data (hb and fs)
            load([path '/FF_' num2str(pair) '_MES_Probe1.mat']);
            hb1(:,1:totCh/2)=hbr; %first 22 channels right hem
            load([path '/FF_' num2str(pair) '_MES_Probe3.mat']);
            hb2(:,1:totCh/2)=hbr;
            
            %left hemisphere both participants
            %get data (hb and fs)
            load([path '/FF_' num2str(pair) '_MES_Probe2.mat']); 
            hb1(:,totCh/2+1:totCh)=hbr; %next 22 channels left hem
            load([path '/FF_' num2str(pair) '_MES_Probe4.mat']); 
            hb2(:,totCh/2+1:totCh)=hbr;
        else
            error('chromophore not specified.');
        end
        
        %now calculate sliding window coherence for every channel
        for channel=1:length(channellist)%Nch
            
            ch=channellist(channel);
            
            %get signal of interest = soi
            if ~isnan(hb1(1, ch)) && ~isnan(hb2(1, ch)) % check if this 
                %channel was not rejected in both subjects during 
                %preprocessing
                sigPart1 = [t, hb1(:,ch)];
                sigPart2 = [t, hb2(:,ch)];
                                
            else
                warning(['Houston, we are having a problem in channel ' ...
                    num2str(ch) 'for dyad ' num2str(pair) ...
                    '. Will not use this channel.']);
                continue; %discard this channel for analysis, try next one.
            end
            
            
            for w=1:nWinTask
               startpoint = task_indices(1)+woffset; %in samples
               part = startpoint+(w-1)*woffset+(w-1)*windowsize: ...
                   startpoint+(w-1)*woffset+w*windowsize-1;                      
                %calc coherence
                c_task_win(w, pair, channel) = ...
                    calculateCoherence(sigPart1, sigPart2, fs, poi, ...
                    coiidx, part);
            end 

            
      % %~~~  for comparison: complete timecourse, online fashion
      
           part = rest_indices; %complete rest block
           c_rest_tot(pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
           
           part = task_indices; %complete task block
           c_task_tot(pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);      
            

    % %~~~  for comparison: complete timecourse, offline analysis
    % 
                %calc coherence
                [pRsq,~, period,coi] = wcoherence(sigPart1(:,2), ...
                    sigPart2(:,2),seconds(1/fs),'PeriodLimits', ...
                    [seconds(6) seconds(14)],'VoicesPerOctave', 14);

                %NaN everything outside cone of influence = coi
                for j=1:1:length(coi)
                    %outer bracket, but for real time need to NaN 
                    %every window at borders. We'll do later.
                    pRsq(period >= coi(j), j) = NaN;
                end
                %reduced pRsq to avoid mayer wave freq
                redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
                %average over periods to get soi
                soi = mean(redpRsq,1,'omitnan'); %soi=signal of interest

                soi_task1 = soi(task_indices);
                soi_rest1 = soi(rest_indices);

                c_rest_tot_offmat(pair, channel ) = mean( ...
                    soi_rest1(coiidx:length(soi_rest1)-coiidx),'omitnan'); 
                c_task_tot_offmat(pair, channel) = mean( ...
                    soi_task1(coiidx:length(soi_task1)-coiidx),'omitnan');

                   
% calculate difference and fill into matrices for later stats   

                counter=0;
                
                for j=1:nWinTask
                    
                    counter=counter+1;
                    
                    diff_task_win(counter, pair, channel) = ...
                        c_task_win(j,pair,channel) - ...
                        c_rest_tot(pair,channel);

                end  
                               
                    diff_task_tot(1,pair,channel) = ...
                        c_task_tot(pair,channel) - ...
                        c_rest_tot(pair,channel);
                    
                
                    diff_task_off(1,pair,channel) = ...
                        c_task_tot_offmat(pair,channel) - ...
                        c_rest_tot_offmat(pair,channel);
                
        
        end %end channelloop
        
    end%end pair loop
    
     disp([num2str(ch) num2str(win) num2str(pair)]);
    % hierarchical binomial model; see function code below.
    teststats_glme(win,:,:) = calculateGlmeForChannels( ...
        diff_task_win, nCh, nPairs, nWinTask);
    all_teststats_glme(win,:,:) = teststats_glme(win,:,:);
    
    
end %end window loop

% %total data, online fashion

teststats_glme_tot = calculateGlmeForChannels( ...
    diff_task_tot, nCh, nPairs, 1); % nWin is one


% %total data, offline fashion

teststats_glme_off = calculateGlmeForChannels( ...
    diff_task_off, nCh, nPairs, 1); % nWin is one

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
    save('workspacemusic_conFB_glme_FDR_TDDR_bPCA_HbO');
    elseif chromophore == 2
    save('workspacemusic_conFB_glme_FDR_TDDR_bPCA_HbR');        
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
end
    er =errorbar(x',bar_data,error_data,'k','linestyle','none', 'LineWidth', 2);
    
    hold off
    
    yline(0,'LineWidth', 2);     
    
    axes=gca;
    axes.FontSize=18;
    axes.FontWeight='b';
    axes.LineWidth=2.5;
    ylim([0 0.15]);
    legend('all channels', 'offline signifcant channels', 'SEM')
    
    
    if chromophore==1
        %title({'difference in WTC, errorbars SEM' , ['for HbO, channel ' ...
        title(['HbO, grandmean over all channels '], ...
            'FontSize', 18, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/Grandmean_HbO_' num2str(ch)], 'svg');
        else
            %do nothing
        end        
    elseif chromophore==2       
        title(['HbR, grandmean over all channels '], ...
            'FontSize', 18, 'FontWeight', 'b');
        if savefig
            saveas(gcf, ['./figures/Grandmean_HbR_' num2str(ch)], 'svg');
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
    end
    
    hold off
    if chromophore==1
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



%% all channels figure


fig=figure;
%fig.Position(3:4)=[1500,2000];
tlt=tiledlayout(10,4,'TileSpacing', 'compact');

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
    if channel < nCh
    X=categorical(1:length(windows));
    elseif channel == nCh
    X=categorical(windows);
    X=reordercats(X,windows); %apparently Matlab is reordering vector when 
    %transforming to categorical... so we need to change back.  
    end
    
    b=bar(X,res_fig_mean(:,channel),'LineWidth', 2.5); hold on
    er = errorbar(X,res_fig_mean(:,channel),res_fig_std(:,channel,2), ...
        res_fig_std(:,channel,1),'LineWidth', 2.5);  
    b.FaceColor = [0.137/0.255, 0.219/0.255, 1];
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    yline(0,'LineWidth', 2.5); 
    
    for i=1:size(all_teststats_glme,1)   
        if ~isnan(all_teststats_glme_FDR(i,channel,9)) 
            text(X(i), er.YData(i)+0.06, '\ast\ast', 'FontSize',22, 'FontWeight', 'b', 'HorizontalAlignment','center');
        elseif ~isnan(all_teststats_glme_FDR(i,channel,8)) 
            text(X(i), er.YData(i)+0.06, '\ast', 'FontSize',22, 'FontWeight', 'b', 'HorizontalAlignment','center');
        else
            %no star
        end
        
    end
    
    axes=gca;
    axes.FontSize=15;%25;
    axes.FontWeight='b';
    axes.LineWidth=3;
    yticks([0 0.1 0.2]);
    ylim([0 0.22]);
    
    hold off

        title(['channel ' ...
            num2str(ch)], 'FontSize', 15, 'FontWeight', 'b');    
    
end


if chromophore==1
    title(tlt,'\Delta WTC per channel and window size for the music data, HbO', 'FontSize', 20, 'FontWeight', 'b');
    if savefig
        saveas(gcf, ['./figures/contFB_allchan_HbO_' num2str(ch)], 'svg');
    else
        %do nothing
    end
elseif chromophore==2
    title(tlt,'\Delta WTC per channel and window size for the music data, HbR', 'FontSize', 20, 'FontWeight', 'b');
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
    rows={'estimate'; 'SE'; 'tStat'; 'DF'; 'pValues'; 'Lower'; 'Upper'; 'q<0.05'};
    mytable=array2table(squeeze(all_teststats_glme_FDR(:,channel,:))', ...
        'VariableNames', windows, 'RowNames', rows);
    
    if chromophore==1
    writetable(mytable,['./tables/music_glme_FDR' num2str(ch) 'HbO_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/music_glme_FDR' num2str(ch) 'HbR_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)   
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

    c_value = mean(soi(coiidx:length(soi) - coiidx), 'omitnan');
end


% mixed model
function teststats_glme = calculateGlmeForChannels(diff_task_win, nCh,... 
    nPairs, nWin)

    teststats_glme = NaN(nCh, 7);

    for channel = 1:nCh
        disp(channel);
        % Extract y values from VZ_tot
       % y = squeeze(VZ_mat(1, :, channel))*nWin;
        y = squeeze(diff_task_win(:, :, channel));
        
        mat = NaN(nPairs * nWin, 2);

        % Populate the result matrix & do Fisher transformation
        for a = 1:nPairs
            for b = 1:nWin
                
                 mat((a - 1) * nWin + b , 1) = atanh(y(b, a));
                 
                 mat((a - 1) *nWin + b, 2) = a;
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
