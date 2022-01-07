%% time windows with music data
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
% compute complete rest block, take blocks of task phase
% taking hbo, as in paper


path=['/Users/kathrin/data/Music_fNIRS/data'];
pairs={'01','02','03','04','05','06', '07','08','09','10','11','12'};
Nch=44; %bad channels will have zeros
phrase=6.61;%sec
%T= 36; %how many windows?  36=about(6600-1800)/fs/(2*phrase) for complete block
%T=32; %31=about(6600-1800-600)/fs/(2*phrase) for block minus learning init
%winsize=round(2*phrase);%roughly 2*phrase=13 sec (multiples of phrase?)
offset=600;
fs=10; %sampling rate 10 Hz
task_indices=1800+600:6600; %signal of interest for task block
%windowsize=300;
rest_indices=1:1800;


windowsizelist=[200 400 600 800 1000 1200];


meanCohen10time=NaN(length(task_indices),length(windowsizelist));
sdCohen10time=NaN(length(task_indices),length(windowsizelist));
meanCohen14time=NaN(length(task_indices),length(windowsizelist));
sdCohen14time=NaN(length(task_indices),length(windowsizelist));
meanCohen31time=NaN(length(task_indices),length(windowsizelist));
sdCohen31time=NaN(length(task_indices),length(windowsizelist));
meanCohen40time=NaN(length(task_indices),length(windowsizelist));
sdCohen40time=NaN(length(task_indices),length(windowsizelist));



for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    defacto_length_results=length(task_indices)-windowsize;
    defacto_length_rest=length(rest_indices)-windowsize;
    
    %coherence
    c_results_rest=NaN(Nch,length(pairs));
    c_results_rest_win=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_task=NaN(Nch,length(pairs),defacto_length_results);
    c_delta=NaN(Nch,length(pairs),defacto_length_results); %matrix with all delta stats per window length
    
    %    c_deltad=zeros(Nch,length(pairs),d,Nwd)*NaN;%for the half trials
    
    
    
    
    %get exact parameters for period and filter
    poi=[6.6 14.02]; %period of interest
    
    load([path '/FF_3_MES_Probe1.mat']); %get data of exemplary pair
    hbo1test=hbo;
    load([path '/FF_3_MES_Probe3.mat']); %get data of exemplary pair
    hbo2test=hbo;
    
    pnoi = zeros(2,1);
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hbo1test(:,14)];
    sigPart2 = [t, hbo2test(:,14)];
    [~,period,~,coi,~] = wtc(sigPart1, sigPart2, 'mcc', 0);
    pnoi(1) = find(period > poi(1), 1, 'first');
    pnoi(2) = find(period < poi(2), 1, 'last');
    
    poi_calc(1)=period(pnoi(1));
    poi_calc(2)=period(pnoi(2));
    
    
    %find cutoff (coi) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    
    
    %    wtc(sigPart1, sigPart2, 'mcc', 0); %have a look at the data
    
    %now real stuff
    for i=1:length(pairs)
        
        %load in participants
        hbo1=zeros(size(hbo1test,1),Nch);
        hbr1=hbo1;
        hbo2=zeros(size(hbo2test,1),Nch);
        hbr2=hbo2;
        %for a start only part learning
        %right hem both participants
        load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
        hbo1(:,1:Nch/2)=hbo; %first 22 channels right hem
        hbr1(:,1:Nch/2)=hbr;
        load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
        hhbo2(:,1:Nch/2)=hbo;
        hbr2(:,1:Nch/2)=hbr;
        
        %left hem both participants
        load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
        hbo1(:,Nch/2+1:Nch)=hbo; %next 22 channels left hem
        hbr1(:,Nch/2+1:Nch)=hbr;
        load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
        hhbo2(:,Nch/2+1:Nch)=hbo;
        hbr2(:,Nch/2+1:Nch)=hbr;
        
        
        parfor ch=1:Nch %for now do this for every channel separately
            
            %get signal of interest
            if ~isnan(hbo1(1, ch)) && ~isnan(hbo2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hbr1(:,ch)];
                sigPart2 = [t, hbr2(:,ch)];
                
                %calc coherence
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % Rsq=r square - measure for coherence
                
                %NaN everything outside coi
                for j=1:1:length(coi)
                    Rsq{ch}(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders
                end
                
                %average over periods to get soi
                soi=mean(Rsq{ch}(pnoi(1):pnoi(2),:),1,'omitnan'); %soi=signal of interest, avaraged over period if interest
            end
            
            %now do window calc
            %calculate rest block completely
            c_results_rest(ch,i)=mean(soi(1:1800),'omitnan');
            
            
            %and different window lengths separately
            
            soi_task=soi(task_indices); %signal of interest for task block
            soi_rest=soi(rest_indices);
            
            c_task_win=mean_wtc_sliding_window(soi_task, cutoff, windowsize);
            c_rest_win=mean_wtc_sliding_window(soi_rest, cutoff, windowsize);
            
            c_results_task(ch,i,:)=c_task_win;
            c_results_rest_win(ch,i,:)=c_rest_win;
            
            %create delta as measure
            if ~isfinite(c_task_win)
                warning(['pair ' num2str(i) 'channel ' num2str(ch) 'needs to be checked.']);
                c_delta(ch,i,:)=NaN(defacto_length_results,1);
            else
                c_delta(ch,i,:)=c_results_task(ch,i,:)-c_results_rest(ch,i);
            end
            
            
            
        end %end of loop over channels
        
        
    end %end of loop over pairs
    
    
    
    %% group stats
    
    %check for infs
    culprits_task=find(~isfinite(c_results_task));
    
    c_delta_noshit=c_delta;
    c_delta_noshit(culprits_task)=NaN;
    
    c_task_noshit=c_results_task;
    c_task_noshit(culprits_task)=NaN;
    
    c_delta_clean=NaN(size(c_delta_noshit));
    %c_corr_clean=NaN(size(c_task_noshit));
    
    c_corr=NaN(Nch,defacto_length_results); %correlation for effect size correction
    c_results_p=NaN(Nch,defacto_length_results);
    Cohens_d_dyad=NaN(Nch, length(pairs), defacto_length_results);
    %Cohens_d_std=NaN(Nch,defacto_length_results);
    %Cohens_d_corrected=NaN(Nch,defacto_length_results);
    
    
    
    for c=1:Nch
        for l=1:defacto_length_results
            for d=1:length(pairs)
                
                %corr for every window for effect size
                %     c_corr(c,l)=corr(squeeze(c_task_noshit(c,:,l))', squeeze(c_results_rest(c,:))');
                
                Cohens_d_dyad(c,d,l)= c_delta(c,d,l)/std(c_results_rest_win(c,d,:));
                
                %    Cohens_d_std(c,l)=std(c_deltaMat,0,1,'omitnan')/sqrt(size(c_deltaMat,1)); %SEM
                %     Cohens_d_corrected(c,l)=Cohens_d_simple(c,l)*1/sqrt(2*(1-c_corr(c,l)));
                
                %         %wilcoxon signed rank for every window
                %
                %         if isempty(find(isnan(c_delta(c,:,l))==0))  %#ok<*EFIND>
                %             continue
                %         end
                %        [p,h,stats] = signrank(squeeze(c_delta(c,:,l)),0,'tail','right');
                %
                %         c_results_p(c,l)=p;
                
                
                %clean data for better display
                %      c_corr_vec_without_outliers=rmoutliers(squeeze(c_corr(c,:,l)),'mean');
                c_delta_vec_without_outliers=rmoutliers(squeeze(c_delta_noshit(c,:,l)),'mean');
                c_delta_clean(c,1:length(c_delta_vec_without_outliers),l)=c_delta_vec_without_outliers;
                
            end
            
        end
    end
    
    
    % c_disp=squeeze(mean(c_delta_clean,2));
    %
    % task_times=(task_indices-(1800+600)+windowsize)/fs;
    % task_times=task_times(1:end-windowsize);
    windowsize_sec=windowsize/fs;
    
    %significant channels in the main study: 10,14,31,40
    %% von hinten auffüllen
    
    meanCohen10time(end-defacto_length_results+1:end,win)=squeeze(mean(Cohens_d_dyad(10,:,:),2));
    sdCohen10time(end-defacto_length_results+1:end,win)=squeeze(std(Cohens_d_dyad(10,:,:),0,2)/sqrt(size(Cohens_d_dyad,2))); %SEM
    meanCohen14time(end-defacto_length_results+1:end,win)=squeeze(mean(Cohens_d_dyad(14,:,:),2));
    sdCohen14time(end-defacto_length_results+1:end,win)=squeeze(std(Cohens_d_dyad(14,:,:),0,2)/sqrt(size(Cohens_d_dyad,2)));
    meanCohen31time(end-defacto_length_results+1:end,win)=squeeze(mean(Cohens_d_dyad(31,:,:),2));
    sdCohen31time(end-defacto_length_results+1:end,win)=squeeze(std(Cohens_d_dyad(31,:,:),0,2)/sqrt(size(Cohens_d_dyad,2)));
    meanCohen40time(end-defacto_length_results+1:end,win)=squeeze(mean(Cohens_d_dyad(40,:,:),2));
    sdCohen40time(end-defacto_length_results+1:end,win)=squeeze(std(Cohens_d_dyad(40,:,:),0,2)/sqrt(size(Cohens_d_dyad,2)));
    
end  %loop over windowsizes
%%

f=figure;
f.Position(3:4)=[1200 900];
tls=tiledlayout('flow');
tls.TileSpacing='none';
tls.Padding='compact';

  title(tls,['Contrast-to-noise ratio (\equiv Cohen''s d)'], 'FontSize', 20, 'FontWeight', 'bold');
  
for  win=1:length(windowsizelist)
    
ax(win)=nexttile;
%x=1+300:3901+300;
%plot(x,Cohens_d_dyad(10,:),x,Cohens_d_dyad(14,:),x,Cohens_d_dyad(40,:),x,Cohens_d_dyad(31,:), 'LineWidth',2);
shadedErrorBar(1:length(task_indices), meanCohen10time(:,win), sdCohen10time(:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color',[0.8500 0.3250 0.0980]},'patchSaturation', 0.3); %red
hold on
shadedErrorBar(1:length(task_indices),meanCohen14time(:,win), sdCohen14time(:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.9290 0.6940 0.1250]},'patchSaturation', 0.3); %yellow 
shadedErrorBar(1:length(task_indices),meanCohen31time(:,win), sdCohen31time(:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.4660 0.6740 0.1880]},'patchSaturation', 0.3); %green
shadedErrorBar(1:length(task_indices),meanCohen40time(:,win), sdCohen40time(:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0 0.4470 0.7410]},'patchSaturation', 0.3); %blue
hold off 
yline(0);

%  ax.XTick=1:600:4000;
%  strval=mat2str(windowsize/fs:60:4000/fs);
%  strval=erase(strval,'[');
%  strval=erase(strval,']');
%  strval=split(strval,' ');
%  ax.XTickLabel=strval;
 ax(win).FontSize = 14;
 %set(Im, 'XData',task_times);

legend('channel 10', 'channel 14', 'channel 40', 'channel 31');
 xlabel('time (samples)');
 ylabel('Contrast-to-noise ratio');
  title(ax(win), ['windowsize=' num2str(windowsizelist(win)/10) ' sec']);
  
end
linkaxes([ax(1),ax(2), ax(3), ax(4), ax(5)],'xy');

%% make window histogram

meanmeanCohen10=mean(meanCohen10time,1,'omitnan');
sdmeanCohen10=mean(sdCohen10time,1,'omitnan');%std(meanCohen10time,0,1,'omitnan')/sqrt(size(meanCohen10time,1));
meanmeanCohen14=mean(meanCohen14time,1,'omitnan');
sdmeanCohen14=mean(sdCohen14time,1,'omitnan');%std(meanCohen14time,0,1,'omitnan')/sqrt(size(meanCohen14time,1));
meanmeanCohen31=mean(meanCohen31time,1,'omitnan');
sdmeanCohen31=mean(sdCohen31time,1,'omitnan');%std(meanCohen31time,0,1,'omitnan')/sqrt(size(meanCohen31time,1));
meanmeanCohen40=mean(meanCohen40time,1,'omitnan');
sdmeanCohen40=mean(sdCohen40time,1,'omitnan');%std(meanCohen40time,0,1,'omitnan')/sqrt(size(meanCohen40time,1));

f=figure;
f.Position(3:4)=[800 600];
tls=tiledlayout('flow');
txt=title(tls,'Contrast-to-noise-ratio (CNR) for different integration window lengths');
txt.FontSize=18;
txt.FontWeight='bold';

windows=windowsizelist/10;

ax1=nexttile;
bar(windows, meanmeanCohen10, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen10,sdmeanCohen10,sdmeanCohen10,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax1.FontSize=13;
ax1.YLim=[-1 10];
ylabel('CNR');
xlabel('window length (sec)');
title(ax1, 'channel 10');

ax2=nexttile;
bar(windows, meanmeanCohen14, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen14,sdmeanCohen14,sdmeanCohen14,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax2.FontSize=13;
ylabel('CNR');
xlabel('window length (sec)');
title(ax2, 'channel 14');

ax3=nexttile;
bar(windows, meanmeanCohen31, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen31,sdmeanCohen31,sdmeanCohen31,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax3.FontSize=13;
ylabel('CNR');
xlabel('window length (sec)');
title(ax3, 'channel 31');

ax4=nexttile;
bar(windows, meanmeanCohen40, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen40,sdmeanCohen40,sdmeanCohen40,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax4.FontSize=13;
ylabel('CNR');
xlabel('window length (sec)');
title(ax4, 'channel 40');

%linkaxes([ax1,ax2, ax3, ax4],'y');
linkaxes([ax2, ax3, ax4],'y');
