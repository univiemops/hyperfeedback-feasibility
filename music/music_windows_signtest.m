%% coherence calculation and sign test for the music data set.
% This script calculates the wtc coherence within a certain frequency band 
%(period of interest) for different task windows.
% It then calculates the difference in coherence between the task
% window/unit of interest and the complete rest block. Statistics are
% assessed via a sign test, indicating if there are significantly more
% windows with higher coherence for task compared to the baseline
% condition. The ratios and numbers of correct vs. incorrect
% classifications are plotted for illustration. This would be ratio of
% correct vs. incorrect feedback in a real experiment.
% 
% KK May 2023



path='YourPathToData/ppData'; %use your data path; Unix System used, adapt accordingly

%experiment specifics
pairs={'01','02','03','04','05','06', '07','08','09','10','11','12'};
nPairs=length(pairs);
totCh=44; %Number of channels; bad channels will have zeros
fs=10; %sampling rate 10 Hz
task_indices=1800+600:6600; %signal of interest for task block
rest_indices=1:1800; %experiment starts with rest

poi=[6 14.0]; %poi = period of interest

%windowlengths to be examined in samples
windowsizelist=[500 600 700 800 900 1000]; %in samples
nTestWindows=length(windowsizelist);


woffset=round(8*fs); %in samples, 8 sec offset between windows


%if only significant channels wanted:
channellist=[10,14,31,40];
nCh=length(channellist);

%initialize variables for results of sign test (coherence_task > coherence_rest) and resulting signs (VZ)
results=NaN(4,nTestWindows,nCh); %(p,zval,sign, total);
results_tot=NaN(4,nCh);
results_off=NaN(4,nCh);
VZ_win=NaN(3, nTestWindows, nPairs, nCh); %3: positive sign, negative + bindings
VZ_tot=NaN(3,nPairs,nCh);
VZ_off=NaN(3,nPairs,nCh);

%chromophore switch
chromophore=2; %1= HbO, 2=HbR
%in paper HbO is used

    
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
[~,~, periodtry,coitry]=wcoherence(sigPart1(:,2), sigPart2(:,2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);

coiidx=find(coitry>max(periodtry),1,'first'); %in samples

%%
for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    %number of windows
    nWinRest=floor(length(rest_indices)/(windowsize+woffset)); %number of possible windows per dyad 
    nWinTask=floor(length(task_indices)/(windowsize+woffset));
    
    %resulting coherence values
    c_rest_win=NaN(nWinRest, nPairs, nCh);
    c_task_win=NaN(nWinTask, nPairs, nCh);

    c_rest_tot=NaN(nPairs, nCh);
    c_task_tot=NaN(nPairs, nCh);
    
    c_rest_tot_offmat=NaN(nPairs, nCh);
    c_task_tot_offmat=NaN(nPairs, nCh);
    
    %matrices for the differences between task and rest coherence
    diff_task_win=NaN(nWinTask, nPairs, nCh);
    diff_task_tot=NaN(nPairs,nCh);
    diff_task_off=NaN(nPairs,nCh);
    
    
    
        %load in pairs, use part learning condition
    for i=1:length(pairs)
        
        %load in participants
        
        hb1=zeros(size(hb1test,1),totCh);
        hb2=zeros(size(hb2test,1),totCh);
        
        
        if chromophore==1
            
            %right hemisphere both participants
            load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
            hb1(:,1:totCh/2)=hbo; %first 22 channels right hem
            load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
            hb2(:,1:totCh/2)=hbo;
            
            %left hemisphere both participants
            load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
            hb1(:,totCh/2+1:totCh)=hbo; %next 22 channels left hem
            load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
            hb2(:,totCh/2+1:totCh)=hbo;
            
        elseif chromophore==2
            
            %right hemisphere both participants
            load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
            hb1(:,1:totCh/2)=hbr; %first 22 channels right hem
            load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
            hb2(:,1:totCh/2)=hbr;
            
            %left hemisphere both participants
            load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
            hb1(:,totCh/2+1:totCh)=hbr; %next 22 channels left hem
            load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
            hb2(:,totCh/2+1:totCh)=hbr;
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
                
                
             % %~~~ calc coherence for win
                
%                    for w=1:nWinRest
%                        startpoint=1+woffset; %in samples
%                        part=startpoint+(w-1)*woffset+(w-1)*windowsize:startpoint+(w-1)*woffset+w*windowsize-1;
%                        tsampvec=startpoint:startpoint+windowsize;
%                        %calc coherence
%                        [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(0.1),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
%                        for j=1:1:length(coi)
%                            pRsq(period >= coi(j), j) = NaN;
%                        end
%                        %get signal of interest
%                        soi=squeeze(mean(pRsq,1,'omitnan'));
%                        soi=soi(part);
%                        c_rest_win(w,i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
%                    end                

                   for w=1:nWinTask
                       startpoint=1800+200+woffset; %in samples
                       part=startpoint+(w-1)*woffset+(w-1)*windowsize:startpoint+(w-1)*woffset+w*windowsize-1;                      
                       %calc coherence
                       [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(0.1),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
                       for j=1:1:length(coi)
                           pRsq(period >= coi(j), j) = NaN;
                       end
                       %get signal of interest
                       soi=squeeze(mean(pRsq,1,'omitnan'));
                       soi=soi(part);
                       c_task_win(w,i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
                   end                    

            else
                warning(['Houston, we are having a problem in channel ' num2str(ch) 'for dyad ' num2str(i) '. Will not use this channel.']);
                continue; %discard this channel for analysis, try next channel.
            end
            
            

            
      % %~~~  for comparison: complete timecourse, online fashion
           part=rest_indices; %complete rest block
           [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(0.1),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
           for j=1:1:length(coi)
               pRsq(period >= coi(j), j) = NaN;
           end
           %get signal of interest
           soi=squeeze(mean(pRsq,1,'omitnan'));
           soi=soi(part);
           c_rest_tot(i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
           
           part=task_indices; %complete task block
           [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(0.1),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
           for j=1:1:length(coi)
               pRsq(period >= coi(j), j) = NaN;
           end
           %get signal of interest
           soi=squeeze(mean(pRsq,1,'omitnan'));
           soi=soi(part);
           c_task_tot(i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');        
            

    % %~~~  for comparison: complete timecourse, offline analysis
    % 
                    %calc coherence
                    [pRsq,~, period,coi]=wcoherence(sigPart1(:,2), sigPart2(:,2),seconds(1/fs),'PeriodLimits', [seconds(6) seconds(14)],'VoicesPerOctave', 14);
    
                    %NaN everything outside cone of influence = coi
                    for j=1:1:length(coi)
                        pRsq(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders. We'll do later.
                    end
    
                    %average over periods to get soi
                    soi=mean(pRsq,1,'omitnan'); %soi=signal of interest, avaraged over period if interest
                    
                    soi_task1=soi(task_indices);
                    soi_rest1=soi(rest_indices);
    
                    c_rest_tot_offmat(i, channel)=mean(soi_rest1(coiidx:length(soi_rest1)-coiidx),'omitnan'); %could also have a look at whole block
                    c_task_tot_offmat(i, channel)=mean(soi_task1(coiidx:length(soi_task1)-coiidx),'omitnan');

            
       
        % ~~~~ count signs

                pos=0;
                neg=0;
                bind=0;
                
                counter=0;
                for j=1:nWinTask
                    counter=counter+1;
                    diff_task_win(counter, i, channel)= c_task_win(j,i,channel)-c_rest_tot(i,channel);
                    if c_task_win(j,i,channel) > c_rest_tot(i,channel)
                        pos=pos+1;
                    elseif c_task_win(j,i,channel)< c_rest_tot(i,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                end
                           
                
                VZ_win(1, win, i, channel)=pos/nWinTask;
                VZ_win(2, win, i, channel)=neg/nWinTask;
                VZ_win(3, win, i, channel)=bind/nWinTask;
        
               %now do the same for the total data online fashion
               
                pos=0;
                neg=0;
                bind=0;
               
                    diff_task_tot(i,channel)=c_task_tot(i,channel)-c_rest_tot(i,channel);
                    if c_task_tot(i,channel)>c_rest_tot(i,channel)
                        pos=pos+1;
                    elseif c_task_tot(i,channel)<c_rest_tot(i,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_tot(1, i, channel)=pos;
                VZ_tot(2, i, channel)=neg;
                VZ_tot(3, i, channel)=bind;
                
            
                
              %and now the same for the offline data  
                pos=0;
                neg=0;
                bind=0;
                
                    diff_task_off(i,channel)=c_task_tot_offmat(i,channel)-c_rest_tot_offmat(i,channel);
                    if c_task_tot_offmat(i,channel)>c_rest_tot_offmat(i,channel)
                        pos=pos+1;
                    elseif c_task_tot_offmat(i,channel)<c_rest_tot_offmat(i,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_off(1, i, channel)=pos;
                VZ_off(2, i, channel)=neg;
                VZ_off(3, i, channel)=bind;  
                
                
                
        
        end %end channelloop
        
    end %end pair loop
    

    %matlab sign test

  for channel=1:nCh

    diffvec=reshape(diff_task_win,[nWinTask*nPairs,nCh]);
    [p,~,stats]= signtest(diffvec(:,channel));
    results(1,win,channel)=p; %(p,zval,sign, total);
    results(2,win,channel)=stats.zval;
    results(3,win,channel)=stats.sign; %number of positive signs
    results(4,win,channel)=nWinTask*nPairs; %total number of signs
    
  end    


end %end window loop

%total data, online fashion
  for channel=1:nCh
    diffvec=reshape(diff_task_tot,[nPairs,nCh]);
    [p,~,stats]= signtest(diffvec(:,channel));
    results_tot(1,channel)=p; %(p,zval,sign, total);
    results_tot(2,channel)=stats.zval;
    results_tot(3,channel)=stats.sign; %number of positive signs
    results_tot(4,channel)=nPairs; %total number of signs    
  end
  
%total data, offline fashion  
  for channel=1:nCh
    diffvec=reshape(diff_task_off,[nPairs,nCh]);
    [p,~,stats]= signtest(diffvec(:,channel));
    results_off(1,channel)=p; %(p,zval,sign, total);
    results_off(2,channel)=stats.zval;
    results_off(3,channel)=stats.sign; %number of positive signs
    results_off(4,channel)=nPairs; %total number of signs
  end
  
  ges_results=results;
  ges_results(:,7,:)=results_tot;
  ges_results(:,8,:)=results_off;
  
%% make total histogram for each channel - ratio of correct vs. incorrectly classified windows per window size

ratios=zeros(size(ges_results,2),nCh);

for channel=1:nCh
     
for i=1:size(ges_results,2)    
ratios(i,channel)=ges_results(3,i,channel)/ges_results(4,i,channel);
        if i<=nTestWindows
            windows(i)={[num2str(round(windowsizelist(i)/fs)) ' sec']};
        elseif i==nTestWindows+1
            windows(i)={'complete'};
        elseif i==nTestWindows+2
            windows(i)={'offline'};
        end
end
figure;
X=categorical(windows);
X=reordercats(X,windows); %apparently Matlab is reordering vector when transforming to categorical... need to change back to correct order.
bar(X,ratios(:,channel),'LineWidth', 4);
axes=gca;
axes.FontSize=25;
axes.FontWeight='b';
axes.LineWidth=4;
axes.YTick=[0.25 0.5 0.75 1];
axes.TickLength = [0.05 0.035];
ylim([0 1.05]);
%xlabel('window size');
if chromophore==1
    title({'Ratio of correctly classified windows' , ['for HbO, channel ' num2str(channellist(channel))]}, 'FontSize', 25, 'FontWeight', 'b');
 %   subtitle(['channel ' num2str(channellist(channel))], 'FontSize', 25, 'FontWeight', 'b');
elseif chromophore==2
    title({'Ratio of correctly classified windows', ['for HbR, channel ' num2str(channellist(channel))]}, 'FontSize', 25, 'FontWeight', 'b');
 %   subtitle(['channel ' num2str(channellist(channel))], 'FontSize', 25, 'FontWeight', 'b');
else 
    warning('Problem with chromophore identity.');
end
  
end

%% write out data

for channel=1:nCh
    %     windows=num2cell(round(windowsizelist/fs));
    ch=channellist(channel);
    for i=1:nTestWindows+2
        if i==nTestWindows+1
            windows(i)={'complete block'};
        elseif i==nTestWindows+2
            windows(i)={'offline'};
        else
            windows(i)={['windowsize ' num2str(round(windowsizelist(i)/fs)) ' sec']};
        end
    end
    rows={'p'; 'z'; 'positive signs'; 'total signs'};
    mytable=array2table(squeeze(ges_results(:,:,channel)), 'VariableNames', windows, 'RowNames', rows);
    
    if chromophore==1
    writetable(mytable,['./tables/music_totsigntest23_channel' num2str(ch) 'HbO.csv'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/music_totsigntest23_channel' num2str(ch) 'HbR.csv'],'WriteRowNames',true)    
    else
        warning('Problem with chromophore identity.');
    end
end


%% bar plot VZ - bonus plot indicating how many pairs have which ratio of correct vs. incorrect windows

fig=figure;
%fig.Position(3:4)=[8000,600];
fig.Position(3:4)=[950,1500];

%mytiles=tiledlayout(nCh, length(windowsizelist)+2, 'TileSpacing', 'compact');
mytiles=tiledlayout(nCh*2, (length(windowsizelist)+2)/2, 'TileSpacing', 'compact');
nTiles=length(windowsizelist)*nCh;
axes=cell(nTiles,1);
tile=1;

for channel=1:nCh
for win=1:length(windowsizelist)+2

    
if win==length(windowsizelist)+1 
    vec=squeeze(VZ_tot(1,:,channel)');  
elseif win==length(windowsizelist)+2
    vec=squeeze(VZ_off(1,:,channel)');  
else
    vec=squeeze(VZ_win(1,win,:,channel));
end

[npair,winratio]=groupcounts(vec);

nexttile
b=bar(winratio, npair, 'FaceColor', 'flat', 'LineWidth', 1.5);
xtips1 = b.XEndPoints;
ytips1 = b.YEndPoints;
labels1 = string(b.YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom', 'FontSize', 14);
axes{tile}=gca;
axes{tile}.FontSize=14;
axes{tile}.LineWidth=1.5;
xticks([0 0.25 0.5 0.75 1]);
xlim([-0.15 1.15])
ylim([0 15]);

if win==length(windowsizelist)+1 
    title(['total, channel ' num2str(channellist(channel))]);
    elseif win==length(windowsizelist)+2
    title(['offline, channel ' num2str(channellist(channel))]);        
else
   title(['win ' num2str(round(windowsizelist(win)/fs)) ' s, channel ' num2str(channellist(channel))]);
end

xlabel('ratio correct feedback');
ylabel('N dyads');
tile=tile+1;
end
end

if chromophore==1
    title(mytiles,'Music: Percentage of correct classification for HbO', 'FontSize', 18, 'FontWeight', 'b');
elseif chromophore==2
    title(mytiles,'Music: Percentage of correct classification for HbR', 'FontSize', 18, 'FontWeight', 'b');
else 
    warning('Problem with chromophore identity.');
end

linkaxes([axes{:}],'xy');