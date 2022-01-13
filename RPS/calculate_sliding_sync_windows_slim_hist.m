       %%  calculate sync values in different windows
    %
    %KK 13.01. now with histogram. To decide: how to pool over rest and/or
    %both task conditions
    %
    %KK 08.01. calc CNR instead of delta wtc
    % note that there is a maximum resolution for calculating CNR
    %since a rest block is only ~65 sec and variability of the window calc
    %in the rest block has to be evaluated. Taking max. ~55 sec (7 trials).
    %taking mean of variability in all three rest blocks for now
    %to do: final figure, group stuff
    %
    %
    %KK 15.12. loop over windowsizes
    %
    %KK 13.10. make nice figures for JGM
    %
    %KK 27.07. improve plot. Need to do: upload debugged get-markers scrpt
    %
    %
    %KK 26.07. create a real sliding window increasing in steps of one
    %sample, debug
    %
    %
    %KK 21.07. make new slimmer script just for wtc windows of double size.
    %removed: restload, windowswitch, assumed avgsub
    %take different starting points
    %
    %KK 30.06. to do: z-score corr already during averaging of rest! check filter for
    %rest!!!!!
    %
    %KK 24.06.21 alternative frequency range: [5 11] - before [6 20]. BUT: would need to
    %recalculate wtc for all rest conditions!
    %
    %KK 08.06.21 including pairs with manual markers now
    %
    %KK 19.05.21 debugging
    %
    %KK 04.05.21 implemented correct coi for every step. rest without
    %correct coi, taken from previous analysis.
    %FDR correction also added to double window.
    %basic bandpassfilter for sliding window corr. also added
    %
    %KK 27.04.21 load in wtc for all rest sessions to get as close as
    %possible to paper result. (load in wtc_rest, do not calculate)
    %to do: tidy up code!
    %
    %KK 19.04.21 adding windowed corr, maybe also need to filter?
    %added mean rest
    %to do: figure out good min window length incorp. coi for each window
    %(prob. about 20-30 sec, about 10-15 sec coi miss), code 'dummy coi miss'
    %also at end of window
    %
    %KK 09.04.21 complete code based on own extracted timestamp mat
    %info: use Grinsted's wavelet toolbox
    %take rest as complete block, tasks in steps of subtask
    %compare rest to following task for different steps


    path=['/Users/kathrin/data/RPS/Kathrin_collab/hmrData'];

    pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};

    conds={'FP','PS', 'PD','C'};

    Nch=16; %number of channels; channel 1-8 DLPFC, channel 9-16 TPJ


    c=1; %using FP for greatest contrast
    d=2; %number of differences of FP-rest per pair 
    T=30;%trials
    f=7.8125; %Hz, sampling rate
    triallength=floor(8*f); %samples per trial, one trial is 8 sec
    RunSamp=8*30*f; %exact number
    RestSamp=509; %65.152 sec rest duration ideally
    %windowlength=triallength*8;%triallength*2; %for now, can change
    
    windowsizelist=[triallength*2 triallength*4 triallength*6];
    channellist=[13,14,15,16];
    
    meanCohentime=NaN(Nch,RunSamp,2,length(windowsizelist));
    sdCohentime=NaN(Nch,RunSamp,2,length(windowsizelist));
    
    meanCohensub=NaN(Nch,length(pairs),2,length(windowsizelist));
    sdCohensub=NaN(Nch,length(pairs),2,length(windowsizelist));    
    
    rmeanCohentime=NaN(Nch,RunSamp, length(windowsizelist));
    rsdCohentime=NaN(Nch, RunSamp, length(windowsizelist));
    
    %% loop over windows
    
    for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    defacto_length_results=RunSamp-windowsize;
    defacto_length_rest=RestSamp-windowsize;
    
    %coherence
 %   c_results_rest=NaN(Nch,length(pairs));
    c_results_rest_win1=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_rest_win2=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_rest_win3=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_task1=NaN(Nch,length(pairs),defacto_length_results);
    c_results_task2=NaN(Nch,length(pairs),defacto_length_results);
    c_delta=NaN(Nch,length(pairs),defacto_length_results*2); %matrix with all delta stats per window length
    

    %for average values
    c_results_rest=zeros(Nch,length(pairs),3)*NaN; 
    mc_results_rest=zeros(Nch,length(pairs))*NaN; %for average over all three rest blocks
    std_avg_rest=NaN(Nch,length(pairs));

    
    c_deltad=zeros(Nch,length(pairs),d,defacto_length_results)*NaN;

    
    %% for coherence - get period of interest poi
    load([path '/Data_' conds{c} '/RPS_01_sub1_' conds{c} '.mat']); %get data of exemplary pair
    hbo1=hbo;
    load([path '/Data_' conds{c} '/RPS_01_sub2_' conds{c} '.mat']); %get data of exemplary pair: hbo, hbr, t (time stamps of measurement)
    hbo2=hbo;
    poi=[5 11];%[6 10];%[5 11]; %[6 20];
    pnoi = zeros(2,1);
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hbo1(:,1)];
    sigPart2 = [t, hbo2(:,1)];
    [~,period,~,coi,~] = wtc(sigPart1, sigPart2, 'mcc', 0);
    pnoi(1) = find(period > poi(1), 1, 'first');
    pnoi(2) = find(period < poi(2), 1, 'last');
    
    poi_calc(1)=period(pnoi(1));
    poi_calc(2)=period(pnoi(2));
    
    
    
    %find cutoff (coi) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    

    %% loop over dyads
    
    for i=1:length(pairs)
        
        
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub1_' conds{c} '.mat']); %get data (hb,t)
        hbo1=hbo;
        hbr1=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub2_' conds{c} '.mat']); %get data (hb,t)
        hbo2=hbo;
        hbr2=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_marker_' conds{c} '.mat']); %get marker matrices
        
        if ~exist('t') || ~exist('fs') ||  ~exist('restMat') ||  ~exist('taskMat') || ~exist('subtaskMat')
            error(['failed to load in markers and/or data for pair ' num2str(pairs{i})]);
        end
        
        %wtc
        %find upper and lower boundaries for period poi
        %get cone of influence coi
        Rsq{Nch} = [];
        Rsq(:) = {NaN(length(period), length(t))};
        
        %% loop over channels of interest
        for chcount=1:length(channellist) %Nch %for now do this for every channel separately
            
            ch=channellist(chcount);
            
            if ~isnan(hbr1(1, ch)) && ~isnan(hbr2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hbr1(:,ch)];
                sigPart2 = [t, hbr2(:,ch)];
                
                %calc coherence
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
                
                %NaN everything outside coi
                for j=1:1:length(coi)
                    Rsq{ch}(period >= coi(j), j) = NaN; %for now use complete calculate, later need coi of substeps
                end
                
                %average over periods to get soi
                soi=mean(Rsq{ch}(pnoi(1):pnoi(2),:),1,'omitnan'); %soi=signal of interest, avaraged over period if interest                
                
            end
            
            
            %now do window calc
            
            %take rest completely for average value
            c_results_rest(ch,i,1)=mean(soi(restMat(1,1):restMat(1,2)),'omitnan');
            c_results_rest(ch,i,2)=mean(soi(restMat(2,1):restMat(2,2)),'omitnan');
            c_results_rest(ch,i,3)=mean(soi(restMat(3,1):restMat(3,2)),'omitnan');
            mc_results_rest(ch,i)=mean(c_results_rest(ch,i,:));
            
            
            soi_task1=soi(subtaskMat(1,1):restMat(2,1)-1);
            soi_task2=soi(subtaskMat(1,2):restMat(3,1)-1);
            soi_rest1=soi(restMat(1,1):restMat(1,2)-1);
            soi_rest2=soi(restMat(2,1):restMat(2,2)-1);
            soi_rest3=soi(restMat(3,1):restMat(3,2)-1);
            
            
         %   % %
%                    mean_wtc_window_vec=NaN(length(soi_task1)-windowsize,1);
% 
%                     for w=1:length(soi_task1)-windowsize
% 
%                         teststat1=soi_task1(w:w+windowsize-1);
% 
%                         for l=1:length(teststat1)
%                             if l<cutoff || l>length(teststat1)-cutoff
%                                 teststat1(l)=NaN;
%                             end
%                         end
%                         mean_wtc_window_vec(w)=mean(teststat1,'omitnan'); 
%                     end
%                     
%                     
            % %
            
            c_task_win1=mean_wtc_sliding_window(soi_task1, cutoff, windowsize);
            c_task_win2=mean_wtc_sliding_window(soi_task2, cutoff, windowsize); 
            
            c_rest_win1=mean_wtc_sliding_window(soi_rest1, cutoff, windowsize); 
            mean_c_rest_win1=mean(c_rest_win1);
            c_rest_win2=mean_wtc_sliding_window(soi_rest2, cutoff, windowsize);   
            mean_c_rest_win2=mean(c_rest_win2);
            c_rest_win3=mean_wtc_sliding_window(soi_rest3, cutoff, windowsize);
            mean_c_rest_win3=mean(c_rest_win3);
            
            %for now only try first rest block
            grandmean_c_rest_win=(mean_c_rest_win1);%+mean_c_rest_win2+ mean_c_rest_win3)/3;  
           
            std_avg_rest(ch,i)=(std(c_rest_win1));%+std(c_rest_win2)+std(c_rest_win3))/3; %since all rest blocks are of equal length, this is the same as mean adjusting each rest block's windows and then calculating total std
                      
            
            c_deltad(ch,i,1,:)=c_task_win1(1:defacto_length_results)-grandmean_c_rest_win;%mc_results_rest(ch,i);  %here are sometimes 2-3 samples missmatch between trigger logging and designed length. Taking designed length to be identical for all sub. For now shouldn't change results but check later!
            c_deltad(ch,i,2,:)=c_task_win2(1:defacto_length_results)-grandmean_c_rest_win;%mc_results_rest(ch,i);
            

        end %end channel loop
    end %end dyad loop
    
%% calc CNR

    Cohens_d_dyad=NaN(Nch, length(pairs),2, defacto_length_results);
    
    for chcount=1:length(channellist)
        ch=channellist(chcount);
        
        for d=1:length(pairs)
            
%            std_rest_avg=(std(c_rest_min1)+std(c_rest_min2)+std(c_rest_min3))/3; %since all rest blocks are of equal length, this is the same as mean adjusting each rest block's windows and then calculating total std
            
            for l=1:defacto_length_results
                for run=1:2
                    
                    %corr for every window for effect size
                    %     c_corr(c,l)=corr(squeeze(c_task_noshit(c,:,l))', squeeze(c_results_rest(c,:))');
                    
                    Cohens_d_dyad(ch,d,run,l)= c_deltad(ch,d,run,l)/std_avg_rest(ch,d);
                    
                    
                end

            end
        end

        for run=1:2

        meanCohentime(ch,end-defacto_length_results+1:end,run,win)=squeeze(mean(Cohens_d_dyad(ch,:,run,:),2));
        sdCohentime(ch,end-defacto_length_results+1:end,run,win)=squeeze(std(Cohens_d_dyad(ch,:,run,:),0,2)/sqrt(size(Cohens_d_dyad,2))); %SEM  
        
        meanCohensub(ch,:,run,win)=squeeze(mean(Cohens_d_dyad(ch,:,run,:),4)); %for trial, shouldn't make a difference
%        sdCohensub=NaN(Nch,length(pairs),2,length(windowsizelist)); 
        
        end
        rmeanCohentime(ch,:,win)=mean(meanCohentime(ch,:,:,win),3);
    %    tmp=squeeze(std(Cohens_d_dyad(ch,:,run,:),0,2));
        rsdCohentime(ch,end-defacto_length_results+1:end,win)=mean(squeeze(std(Cohens_d_dyad(ch,:,run,:),0,2)),3)/sqrt(size(Cohens_d_dyad,2));
    end

    
    
    
    end %end window loop
    
    
    
     %% plot stuff single
     

fig=figure;
fig.Position(3:4)=[1200 900];
tls=tiledlayout('flow');
tls.TileSpacing='none';
tls.Padding='compact';

  title(tls,['Contrast-to-noise ratio (\equiv Cohen''s d)'], 'FontSize', 20, 'FontWeight', 'bold');
  
for  win=1:length(windowsizelist)
    
ax(win)=nexttile;
%x=1+300:3901+300;
%plot(x,Cohens_d_dyad(10,:),x,Cohens_d_dyad(14,:),x,Cohens_d_dyad(40,:),x,Cohens_d_dyad(31,:), 'LineWidth',2);
shadedErrorBar(1:RunSamp, meanCohentime(13,:,1,win), sdCohentime(13,:,1,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color',[0.8500 0.3250 0.0980]},'patchSaturation', 0.3); %red
hold on
shadedErrorBar(1:RunSamp, meanCohentime(14,:,1,win), sdCohentime(14,:,1,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.9290 0.6940 0.1250]},'patchSaturation', 0.3); %yellow 
shadedErrorBar(1:RunSamp, meanCohentime(15,:,1,win), sdCohentime(15,:,1,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.4660 0.6740 0.1880]},'patchSaturation', 0.3); %green
shadedErrorBar(1:RunSamp ,meanCohentime(16,:,1,win), sdCohentime(16,:,1,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0 0.4470 0.7410]},'patchSaturation', 0.3); %blue
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

legend('channel 13', 'channel 14', 'channel 15', 'channel 16');
 xlabel('time (samples)');
 ylabel('Contrast-to-noise ratio');
  title(ax(win), ['windowsize=' num2str(windowsizelist(win)/f) ' sec']);
  
end
%linkaxes([ax(1),ax(2), ax(3), ax(4), ax(5)],'xy');

% %% plot stuff mean
% 
% fig=figure;
% fig.Position(3:4)=[1200 900];
% tls=tiledlayout('flow');
% tls.TileSpacing='none';
% tls.Padding='compact';
% 
%   title(tls,['Contrast-to-noise ratio (\equiv Cohen''s d)'], 'FontSize', 20, 'FontWeight', 'bold');
%   
% for  win=1:length(windowsizelist)
%     
% ax(win)=nexttile;
% %x=1+300:3901+300;
% %plot(x,Cohens_d_dyad(10,:),x,Cohens_d_dyad(14,:),x,Cohens_d_dyad(40,:),x,Cohens_d_dyad(31,:), 'LineWidth',2);
% shadedErrorBar(1:RunSamp, rmeanCohentime(13,:,win), rsdCohentime(13,:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color',[0.8500 0.3250 0.0980]},'patchSaturation', 0.3); %red
% hold on
% shadedErrorBar(1:RunSamp, rmeanCohentime(14,:,win), rsdCohentime(14,:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.9290 0.6940 0.1250]},'patchSaturation', 0.3); %yellow 
% shadedErrorBar(1:RunSamp, rmeanCohentime(15,:,win), rsdCohentime(15,:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0.4660 0.6740 0.1880]},'patchSaturation', 0.3); %green
% shadedErrorBar(1:RunSamp ,rmeanCohentime(16,:,win), rsdCohentime(16,:,win), 'lineProps',{'-', 'LineWidth',1.5, 'Color', [0 0.4470 0.7410]},'patchSaturation', 0.3); %blue
% hold off 
% yline(0);
% 
% %  ax.XTick=1:600:4000;
% %  strval=mat2str(windowsize/fs:60:4000/fs);
% %  strval=erase(strval,'[');
% %  strval=erase(strval,']');
% %  strval=split(strval,' ');
% %  ax.XTickLabel=strval;
%  ax(win).FontSize = 14;
%  %set(Im, 'XData',task_times);
% 
% legend('channel 13', 'channel 14', 'channel 15', 'channel 16');
%  xlabel('time (samples)');
%  ylabel('Contrast-to-noise ratio');
%   title(ax(win), ['windowsize=' num2str(windowsizelist(win)/f) ' sec']);
%   
% end

%% make window histogram

meanmeanCohen13=squeeze(mean(meanCohentime(13,:,1,:),2,'omitnan'));
sdmeanCohen13=squeeze(mean(sdCohentime(13,:,1,:),2,'omitnan'));%std(meanCohen10time,0,1,'omitnan')/sqrt(size(meanCohen10time,1));
meanmeanCohen14=squeeze(mean(meanCohentime(14,:,1,:),2,'omitnan'));
sdmeanCohen14=squeeze(mean(sdCohentime(14,:,1,:),2,'omitnan'));
meanmeanCohen15=squeeze(mean(meanCohentime(15,:,1,:),2,'omitnan'));
sdmeanCohen15=squeeze(mean(sdCohentime(15,:,1,:),2,'omitnan'));
meanmeanCohen16=squeeze(mean(meanCohentime(16,:,1,:),2,'omitnan'));
sdmeanCohen16=squeeze(mean(sdCohentime(16,:,1,:),2,'omitnan'));

fhist=figure;
fhist.Position(3:4)=[800 600];
tls=tiledlayout('flow');
txt=title(tls,'Contrast-to-noise-ratio (CNR) for different integration window lengths');
txt.FontSize=18;
txt.FontWeight='bold';

windows=windowsizelist;

ax1=nexttile;
bar(windows, meanmeanCohen13, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen13,sdmeanCohen13,sdmeanCohen13,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax1.FontSize=13;
ax1.YLim=[-1 10];
ylabel('CNR');
xlabel('window length (sec)');
title(ax1, 'channel 13');

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
bar(windows, meanmeanCohen15, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen15,sdmeanCohen15,sdmeanCohen15,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax3.FontSize=13;
ylabel('CNR');
xlabel('window length (sec)');
title(ax3, 'channel 15');

ax4=nexttile;
bar(windows, meanmeanCohen16, 'LineWidth', 1.5); hold on
er=errorbar(windows,meanmeanCohen16,sdmeanCohen16,sdmeanCohen16,'LineWidth', 1.5);
er.Color=[0 0 0];
er.LineStyle='none';
ax4.FontSize=13;
ylabel('CNR');
xlabel('window length (sec)');
title(ax4, 'channel 16');

%linkaxes([ax1,ax2, ax3, ax4],'y');
linkaxes([ax2, ax3, ax4],'y');



% %% make window histogram - try out other stuff
% 
% meanmeanCohen13=squeeze(mean(meanCohensub(13,:,1,:),2,'omitnan'));
% sdmeanCohen13=squeeze(mean(sdCohensub(13,:,1,:),2,'omitnan'));%std(meanCohen10time,0,1,'omitnan')/sqrt(size(meanCohen10time,1));
% meanmeanCohen14=squeeze(mean(meanCohensub(14,:,1,:),2,'omitnan'));
% sdmeanCohen14=squeeze(mean(sdCohensub(14,:,1,:),2,'omitnan'));
% meanmeanCohen15=squeeze(mean(meanCohensub(15,:,1,:),2,'omitnan'));
% sdmeanCohen15=squeeze(mean(sdCohensub(15,:,1,:),2,'omitnan'));
% meanmeanCohen16=squeeze(mean(meanCohensub(16,:,1,:),2,'omitnan'));
% sdmeanCohen16=squeeze(mean(sdCohensub(16,:,1,:),2,'omitnan'));
% 
% fhist=figure;
% fhist.Position(3:4)=[800 600];
% tls=tiledlayout('flow');
% txt=title(tls,'Contrast-to-noise-ratio (CNR) for different integration window lengths');
% txt.FontSize=18;
% txt.FontWeight='bold';
% 
% windows=windowsizelist;
% 
% ax1=nexttile;
% bar(windows, meanmeanCohen13, 'LineWidth', 1.5); hold on
% er=errorbar(windows,meanmeanCohen13,sdmeanCohen13,sdmeanCohen13,'LineWidth', 1.5);
% er.Color=[0 0 0];
% er.LineStyle='none';
% ax1.FontSize=13;
% ax1.YLim=[-1 10];
% ylabel('CNR');
% xlabel('window length (sec)');
% title(ax1, 'channel 13');
% 
% ax2=nexttile;
% bar(windows, meanmeanCohen14, 'LineWidth', 1.5); hold on
% er=errorbar(windows,meanmeanCohen14,sdmeanCohen14,sdmeanCohen14,'LineWidth', 1.5);
% er.Color=[0 0 0];
% er.LineStyle='none';
% ax2.FontSize=13;
% ylabel('CNR');
% xlabel('window length (sec)');
% title(ax2, 'channel 14');
% 
% ax3=nexttile;
% bar(windows, meanmeanCohen15, 'LineWidth', 1.5); hold on
% er=errorbar(windows,meanmeanCohen15,sdmeanCohen15,sdmeanCohen15,'LineWidth', 1.5);
% er.Color=[0 0 0];
% er.LineStyle='none';
% ax3.FontSize=13;
% ylabel('CNR');
% xlabel('window length (sec)');
% title(ax3, 'channel 15');
% 
% ax4=nexttile;
% bar(windows, meanmeanCohen16, 'LineWidth', 1.5); hold on
% er=errorbar(windows,meanmeanCohen16,sdmeanCohen16,sdmeanCohen16,'LineWidth', 1.5);
% er.Color=[0 0 0];
% er.LineStyle='none';
% ax4.FontSize=13;
% ylabel('CNR');
% xlabel('window length (sec)');
% title(ax4, 'channel 16');
% 
% %linkaxes([ax1,ax2, ax3, ax4],'y');
% linkaxes([ax2, ax3, ax4],'y');

