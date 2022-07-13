%% time windows with music data
%
%12.07. clean script, add switches
%
%11.07. bugfixes - totCohen for correct total effect size for bar plots
%
%25.03. minor edits for figure
%
%15.10.21
%
%KK 17.08. adapting for sliding wtc windows. throwing out corr for now.
%
%KK 29.06., 30.06. adapting from calculate_sync_windows for RPS
%beware: check faulty channels from preprocessing!
% to do: filter for corr
%
%
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
%significant channels in the main study: 10,14,31,40


path=['/Users/kathrin/data/Music_fNIRS/data']; %use your data path

%experiment specifics
pairs={'01','02','03','04','05','06', '07','08','09','10','11','12'};
Nch=44; %Number of channels; bad channels will have zeros
fs=10; %sampling rate 10 Hz
task_indices=1800+600:6600; %signal of interest for task block
rest_indices=1:1800; %experiment starts with rest

poi=[6.6 14.02]; %poi = period of interest

%windowlengths to be examined in samples
windowsizelist=[200 300 400 500];

%if only significant channels wanted:
channellist=[10,14,31,40];

%initialize all kinds of variables
meanCohen10time=NaN(length(task_indices),length(windowsizelist));
sdCohen10time=NaN(length(task_indices),length(windowsizelist));
meanCohen14time=NaN(length(task_indices),length(windowsizelist));
sdCohen14time=NaN(length(task_indices),length(windowsizelist));
meanCohen31time=NaN(length(task_indices),length(windowsizelist));
sdCohen31time=NaN(length(task_indices),length(windowsizelist));
meanCohen40time=NaN(length(task_indices),length(windowsizelist));
sdCohen40time=NaN(length(task_indices),length(windowsizelist));

totmeanCohen10=NaN(length(windowsizelist),1);
totsdCohen10=NaN(length(windowsizelist),1);
totmeanCohen14=NaN(length(windowsizelist),1);
totsdCohen14=NaN(length(windowsizelist),1);
totmeanCohen31=NaN(length(windowsizelist),1);
totsdCohen31=NaN(length(windowsizelist),1);
totmeanCohen40=NaN(length(windowsizelist),1);
totsdCohen40=NaN(length(windowsizelist),1);

%chromophore switch
chromophore=2; %1= HbO, 2=HbR


%loop over windowsizes
for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    defacto_length_results=length(task_indices)-windowsize;
    defacto_length_rest=length(rest_indices)-windowsize;
    
    %initialize coherence variables
    c_results_rest_win=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_task=NaN(Nch,length(pairs),defacto_length_results);
    
    mc_results_rest_win=NaN(Nch,length(pairs));
    mc_results_task=NaN(Nch,length(pairs));
    
    c_delta=NaN(Nch,length(pairs),defacto_length_results); %matrix with all delta stats per window length
    mc_delta=NaN(Nch,length(pairs));
    
    
    
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
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hb1test(:,14)];
    sigPart2 = [t, hb2test(:,14)];
    [~,period,~,coi,~] = wtc(sigPart1, sigPart2, 'mcc', 0);
    pnoi(1) = find(period > poi(1), 1, 'first');
    pnoi(2) = find(period < poi(2), 1, 'last');
    
    poi_calc(1)=period(pnoi(1));
    poi_calc(2)=period(pnoi(2));
    
    
    %find cutoff (coi = cone of influence) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    
    
    
    %now do the real stuff!
    
    %load in pairs, use part learning condition
    for i=1:length(pairs)
        
        %load in participants
        
        hb1=zeros(size(hb1test,1),Nch);
        %       hbr1=zeros(size(hb1test,1),Nch);
        hb2=zeros(size(hb2test,1),Nch);
        %      hbr2=zeros(size(hb2test,1),Nch);
        
        
        if chromophore==1
            
            %right hem both participants
            load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
            hb1(:,1:Nch/2)=hbo; %first 22 channels right hem
            load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
            hb2(:,1:Nch/2)=hbo;
            
            %left hem both participants
            load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
            hb1(:,Nch/2+1:Nch)=hbo; %next 22 channels left hem
            load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
            hb2(:,Nch/2+1:Nch)=hbo;
            
        elseif chromophore==2
            
            %right hem both participants
            load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
            hb1(:,1:Nch/2)=hbr; %first 22 channels right hem
            load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
            hb2(:,1:Nch/2)=hbr;
            
            %left hem both participants
            load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
            hb1(:,Nch/2+1:Nch)=hbr; %next 22 channels left hem
            load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
            hb2(:,Nch/2+1:Nch)=hbr;
        else
            error('chromophore not specified.');
        end
        
        %now calculate sliding window coherence for every channel
        for channel=1:length(channellist)%Nch
            
            ch=channellist(channel);
            
            %get signal of interest = soi
            if ~isnan(hb1(1, ch)) && ~isnan(hb2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hb1(:,ch)];
                sigPart2 = [t, hb2(:,ch)];
                
                %calc coherence
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);    % Rsq = r square - measure for coherence
                
                %NaN everything outside cone of influence = coi
                for j=1:1:length(coi)
                    Rsq{ch}(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders. We'll do later.
                end
                
                %average over periods to get soi
                soi=mean(Rsq{ch}(pnoi(1):pnoi(2),:),1,'omitnan'); %soi=signal of interest, avaraged over period if interest
            else
                warning(['houston, we are having a problem in channel ' num2str(ch) 'for dyad ' num2str(i)]);
            end
            
            
            %now do window calc
            soi_task=soi(task_indices); %soi for task block
            soi_rest=soi(rest_indices); %soi for rest block
            
            %and finally: here's the sliding window
            c_task_win=mean_wtc_sliding_window(soi_task, cutoff, windowsize);
            c_rest_win=mean_wtc_sliding_window(soi_rest, cutoff, windowsize);
            
            %assign to mats
            c_results_task(ch,i,:)=c_task_win; %we'll use this later for sliding window CNR = contrast-to-noise ratio
            mc_results_task(ch,i)=mean(c_results_task(ch,i,:), 'omitnan');
            c_results_rest_win(ch,i,:)=c_rest_win; %not used, only for calculation of mean
            mc_results_rest_win(ch,i)=mean(c_results_rest_win(ch,i,:), 'omitnan');
            
            %create delta as measure
            if any(~isfinite(c_task_win)) || any(~isfinite(c_rest_win))
                warning(['pair ' num2str(i) 'channel ' num2str(ch) 'needs to be checked.']);
                c_delta(ch,i,:)=NaN(defacto_length_results,1);
                mc_delta(ch,i)=NaN;
            else
                c_delta(ch,i,:)=c_results_task(ch,i,:)-mc_results_rest_win(ch,i); %for sliding window CNR plots
                mc_delta(ch,i)=mc_results_task(ch,i)-mc_results_rest_win(ch,i); %for histograms
            end
            %
            
        end %end of loop over channels
        
        
    end %end of loop over pairs
    
    
    
    %% group stats
    
    
    %initialize
    winCohens_d_dyad=NaN(Nch, length(pairs), defacto_length_results); %for sliding windowed CNR
    totCohens_d_dyad=NaN(Nch, length(pairs)); %for histograms
    
    
    for c=1:Nch
        for l=1:defacto_length_results
            for d=1:length(pairs)
                
                winCohens_d_dyad(c,d,l)= c_delta(c,d,l)/std(c_results_rest_win(c,d,:));
                totCohens_d_dyad(c,d)=mc_delta(c,d)/std(c_results_rest_win(c,d,:));
            end
            
        end
    end
    
    
    %% write out data, start from back side
    %significant channels in the main study: 10,14,31,40
    windowsize_sec=windowsize/fs;
    
    %for sliding CNR plots
    meanCohen10time(end-defacto_length_results+1:end,win)=squeeze(mean(winCohens_d_dyad(10,:,:),2));
    sdCohen10time(end-defacto_length_results+1:end,win)=squeeze(std(winCohens_d_dyad(10,:,:),0,2)/sqrt(size(winCohens_d_dyad,2))); %SEM
    meanCohen14time(end-defacto_length_results+1:end,win)=squeeze(mean(winCohens_d_dyad(14,:,:),2));
    sdCohen14time(end-defacto_length_results+1:end,win)=squeeze(std(winCohens_d_dyad(14,:,:),0,2)/sqrt(size(winCohens_d_dyad,2)));
    meanCohen31time(end-defacto_length_results+1:end,win)=squeeze(mean(winCohens_d_dyad(31,:,:),2));
    sdCohen31time(end-defacto_length_results+1:end,win)=squeeze(std(winCohens_d_dyad(31,:,:),0,2)/sqrt(size(winCohens_d_dyad,2)));
    meanCohen40time(end-defacto_length_results+1:end,win)=squeeze(mean(winCohens_d_dyad(40,:,:),2));
    sdCohen40time(end-defacto_length_results+1:end,win)=squeeze(std(winCohens_d_dyad(40,:,:),0,2)/sqrt(size(winCohens_d_dyad,2)));
    
    %total "static" CNR, for histograms and reporting
    totmeanCohen10(win)=squeeze(mean(totCohens_d_dyad(10,:),2));
    totsdCohen10(win)=squeeze(std(totCohens_d_dyad(10,:),0,2)/sqrt(size(totCohens_d_dyad,2))); %SEM
    totmeanCohen14(win)=squeeze(mean(totCohens_d_dyad(14,:),2));
    totsdCohen14(win)=squeeze(std(totCohens_d_dyad(14,:),0,2)/sqrt(size(totCohens_d_dyad,2))); %SEM
    totmeanCohen31(win)=squeeze(mean(totCohens_d_dyad(31,:),2));
    totsdCohen31(win)=squeeze(std(totCohens_d_dyad(31,:),0,2)/sqrt(size(totCohens_d_dyad,2))); %SEM
    totmeanCohen40(win)=squeeze(mean(totCohens_d_dyad(40,:),2));
    totsdCohen40(win)=squeeze(std(totCohens_d_dyad(40,:),0,2)/sqrt(size(totCohens_d_dyad,2))); %SEM
    
end  %loop over windowsizes



%% Plot data, sliding CNR

fig=figure;
fig.Position(3:4)=[1800 900];
tls=tiledlayout('flow');
tls.TileSpacing='none';
tls.Padding='compact';

title(tls,['Contrast-to-noise ratio (\equiv Cohen''s d), parts teaching, r+l IFG/DLPFC'], 'FontSize', 20, 'FontWeight', 'bold');


for  win=1:length(windowsizelist)
    
    ax(win)=nexttile;
    
    shadedErrorBar((1:length(task_indices))/fs, meanCohen10time(:,win), sdCohen10time(:,win), 'lineProps',{'-', 'LineWidth',2, 'Color',[0.8500 0.3250 0.0980]},'patchSaturation', 0.3); %red
    hold on
    shadedErrorBar((1:length(task_indices))/fs,meanCohen14time(:,win), sdCohen14time(:,win), 'lineProps',{'-', 'LineWidth',2, 'Color', [0.9290 0.6940 0.1250]},'patchSaturation', 0.3); %yellow
    shadedErrorBar((1:length(task_indices))/fs,meanCohen31time(:,win), sdCohen31time(:,win), 'lineProps',{'-', 'LineWidth',2, 'Color', [0.4660 0.6740 0.1880]},'patchSaturation', 0.3); %green
    shadedErrorBar((1:length(task_indices))/fs,meanCohen40time(:,win), sdCohen40time(:,win), 'lineProps',{'-', 'LineWidth',2, 'Color', [0 0.4470 0.7410]},'patchSaturation', 0.3); %blue
    hold off
    yline(0, 'LineWidth', 1.5);
    ax(win).LineWidth = 1.5;
    box on;
    
    ax(win).FontSize = 18;
    
    legend('channel 10', 'channel 14', 'channel 40', 'channel 31');
    xlabel('time (sec)');
    ylabel('Contrast-to-noise ratio');
    title(ax(win), ['windowsize=' num2str(windowsizelist(win)/10) ' sec']);
    ax(win).XLim=[0 420];
end
linkaxes([ax(1),ax(2), ax(3), ax(4)],'xy');



%% make window histogram for absolute CNR

fhist=figure;
fhist.Position(3:4)=[800 600];
tls=tiledlayout('flow');
txt=title(tls,'music data: Contrast-to-noise-ratio (CNR) for different integration window lengths');
txt.FontSize=18;
txt.FontWeight='bold';

windows=windowsizelist/10;

ax1=nexttile;
b=bar(windows, totmeanCohen10, 'LineWidth', 1.5); hold on
b.BaseLine.LineWidth=1.5;
er=errorbar(windows,totmeanCohen10,totsdCohen10,totsdCohen10,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax1.FontSize=18;
ax1.LineWidth = 1.5;
box on;
ax1.YLim=[-0.5 3];
ylabel('CNR');
xlabel('window length (sec)');
title(ax1, 'channel 10');

ax2=nexttile;
b=bar(windows, totmeanCohen14, 'LineWidth', 1.5); hold on
b.BaseLine.LineWidth=1.5;
er=errorbar(windows,totmeanCohen14,totsdCohen14,totsdCohen14,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax2.FontSize=18;
ax2.LineWidth = 1.5;
box on;
ylabel('CNR');
xlabel('window length (sec)');
title(ax2, 'channel 14');

ax3=nexttile;
b=bar(windows, totmeanCohen31, 'LineWidth', 1.5); hold on
b.BaseLine.LineWidth=1.5;
er=errorbar(windows,totmeanCohen31,totsdCohen31,totsdCohen31,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax3.FontSize=18;
ax3.LineWidth = 1.5;
box on;
ylabel('CNR');
xlabel('window length (sec)');
title(ax3, 'channel 31');

ax4=nexttile;
b=bar(windows, totmeanCohen40, 'LineWidth', 1.5); hold on
b.BaseLine.LineWidth=1.5;
er=errorbar(windows,totmeanCohen40,totsdCohen40,totsdCohen40,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax4.FontSize=18;
ax4.LineWidth = 1.5;
box on;
ylabel('CNR');
xlabel('window length (sec)');
title(ax4, 'channel 40');

linkaxes([ax1, ax2, ax3, ax4],'xy');
