%% coherence calculation and sign test for RPS data set.
% This script calculates the wavelet tranform coherence within a certain frequency band 
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


path=['yourPathToData/hmrData']; %Unix System; change if on Windows

pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};
nPairs=length(pairs);

totCh=16; %number of channels; channel 1-8 DLPFC, channel 9-16 TPJ

fs=7.8125; %Hz, sampling rate
T=30;%trials
triallength=floor(8*fs); %samples per trial, one trial is 8 sec
RunSamp=8*30*fs; %exact number
RestSamp=509; %65.152 sec rest duration ideally

conds={'FP','PS', 'PD','C'};
co=1; %using FP for greatest contrast

windowsizelist=[round(50*fs) round(60*fs) round(70*fs) round(80*fs) round(90*fs) round(100*fs)]; %in samples

nTestWindows=length(windowsizelist);

woffset=round(8*fs); %in samples, offset between windows

channellist=[13,14,15,16];
nCh=length(channellist);

poi=[6 14];

%chromophore switch
chromophore=2; %1= HbO, 2=HbR; HbR is used in Kayhan et al. (2022)

%initialize variables for results of sign test (coherence_task > coherence_rest) and resulting signs (VZ)
results=NaN(4,nTestWindows,nCh); %(p,zval,sign, total);
results_tot=NaN(4,nCh);
results_off=NaN(4,nCh);
VZ_win=NaN(3, nTestWindows, nPairs, nCh); %3: positive sign, negative + bindings
VZ_tot=NaN(3,nPairs,nCh);
VZ_off=NaN(3,nPairs,nCh);



%get exact parameters for period and filter

%get data of exemplary pair
if chromophore==1
    load([path '/Data_' conds{co} '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test=hbo;
    load([path '/Data_' conds{co} '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test=hbo;
elseif chromophore==2
    load([path '/Data_' conds{co} '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test=hbr;
    load([path '/Data_' conds{co} '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test=hbr;
else
    error('chromophore not specified.');
end

%get exemplary markers
load([path '/Data_' conds{co} '/RPS_01_marker_' conds{co} '.mat']); %get marker matrices
if ~exist('t') || ~exist('fs') ||  ~exist('restMat') ||  ~exist('taskMat') || ~exist('subtaskMat')
    error(['failed to load in markers and/or data for pair1']);
end

task_indices1=subtaskMat(1,1):restMat(2,1)-1;
task_indices2=subtaskMat(1,2):restMat(3,1)-1;
rest_indices1=restMat(1,1):restMat(1,2)-1;
rest_indices2=restMat(2,1):restMat(2,2)-1;
rest_indices3=restMat(3,1):restMat(3,2)-1;


sigPart1 = [t, hb1test(:,14)];
sigPart2 = [t, hb2test(:,14)];


%get cutoff for coi for Matlab wavelet toolbox:
[~,~, periodtry,coitry]=wcoherence(sigPart1(:,2), sigPart2(:,2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);

coiidx=find(coitry>max(periodtry),1,'first'); %in samples



for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    %number of windows
    nWinRest1=floor(length(rest_indices1)/(windowsize+woffset));
    nWinRest2=floor(length(rest_indices2)/(windowsize+woffset));
    nWinRest3=floor(length(rest_indices3)/(windowsize+woffset));
    
    nWinTask1=floor(length(task_indices1)/(windowsize+woffset));
    nWinTask2=floor(length(task_indices2)/(windowsize+woffset));
    
    nWinRest=nWinRest1+nWinRest2+nWinRest3;
    nWinTask=nWinTask1+nWinTask2;
    
    %resulting coherence values
    c_rest_win=NaN(nWinRest, nPairs, nCh);
    c_task_win=NaN(nWinTask, nPairs, nCh);
    
    c_rest_tot=NaN(3,nPairs, nCh);
    c_task_tot=NaN(2,nPairs, nCh);
 
    c_rest_tot_offmat=NaN(3,nPairs, nCh);
    c_task_tot_offmat=NaN(2,nPairs, nCh);
    
    %matrices for the differences between task and rest coherence
    diff_task_win=NaN(nWinTask, nPairs, nCh);
    diff_task_tot=NaN(2,nPairs,nCh);
    diff_task_off=NaN(2,nPairs,nCh);
    


    %load in pairs, use part learning condition
    for i=1:length(pairs)
        
        %load in participants
        
        if chromophore==1
            
            load([path '/Data_' conds{co} '/RPS_' pairs{i} '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1=hbo;
            load([path '/Data_' conds{co} '/RPS_' pairs{i} '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2=hbo;
            
        elseif chromophore==2
            
            load([path '/Data_' conds{co} '/RPS_' pairs{i} '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1=hbr;
            load([path '/Data_' conds{co} '/RPS_' pairs{i} '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2=hbr;
            
        else
            error('chromophore not specified.');
        end
        
        %in case markers differ slightly for this subject, extract them again:
        load([path '/Data_' conds{co} '/RPS_' pairs{i} '_marker_' conds{co} '.mat']); %get marker matrices
        if ~exist('t') || ~exist('fs') ||  ~exist('restMat') ||  ~exist('taskMat') || ~exist('subtaskMat')
            error(['failed to load in markers and/or data for pair ' num2str(pairs{i})]);
        end
        task_indices1=subtaskMat(1,1):restMat(2,1)-1;
        task_indices2=subtaskMat(1,2):restMat(3,1)-1;
        rest_indices1=restMat(1,1):restMat(1,2)-1;
        rest_indices2=restMat(2,1):restMat(2,2)-1;
        rest_indices3=restMat(3,1):restMat(3,2)-1;
        
        
        %~~~~ now calculate coherence for every channel
        for channel=1:length(channellist)%Nch
            
            ch=channellist(channel);
            
            %get signal of interest = soi
            if ~isnan(hb1(1, ch)) && ~isnan(hb2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hb1(:,ch)];
                sigPart2 = [t, hb2(:,ch)];
                
                
                %~~~~ online windowed coherence in windows
                
                w_counter=0;
                for w=1:nWinTask1
                    w_counter=w_counter+1;
                    startpoint=task_indices1(1)+woffset; %in samples
                    part=startpoint+(w-1)*woffset+(w-1)*windowsize:startpoint+(w-1)*woffset+w*windowsize-1;
                    %calc coherence
                    [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
                    for j=1:1:length(coi) %NaN everything outside of coi. just a safety net, windowed application of coi below.
                        pRsq(period >= coi(j), j) = NaN;
                    end
                    %get signal of interest =soi
                    soi=squeeze(mean(pRsq,1,'omitnan'));
                    soi=soi(part);
                    c_task_win(w_counter,i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
                end
                for w=1:nWinTask2
                    w_counter=w_counter+1;
                    startpoint=task_indices2(1)+woffset; %in samples
                    part=startpoint+(w-1)*woffset+(w-1)*windowsize:startpoint+(w-1)*woffset+w*windowsize-1;
                    %calc coherence
                    [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
                    for j=1:1:length(coi)
                        pRsq(period >= coi(j), j) = NaN;
                    end
                    %get signal of interest
                    soi=squeeze(mean(pRsq,1,'omitnan'));
                    soi=soi(part);
                    c_task_win(w_counter,i,channel)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
                end
                
            else
                warning(['houston, we are having a problem in channel ' num2str(ch) 'for dyad ' num2str(i)]);
            end
            
            
            
            %~~~~ online complete blocks
            
            %rest blocks
            part=rest_indices1;
            [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            soi1=squeeze(mean(pRsq,1,'omitnan'));
            soi1=soi1(rest_indices1);
            c_rest_tot(1,i,channel)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
            
            part=rest_indices2;
            [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            soi2=squeeze(mean(pRsq,1,'omitnan'));
            soi2=soi2(rest_indices2);
            c_rest_tot(2,i,channel)=mean(soi2(coiidx:length(soi2)-coiidx),'omitnan');
            
            part=rest_indices3;
            [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            soi3=squeeze(mean(pRsq,1,'omitnan'));
            soi3=soi3(rest_indices3);
            c_rest_tot(3,i,channel)=mean(soi3(coiidx:length(soi3)-coiidx),'omitnan');
            
            
            %task blocks
            part=task_indices1;
            [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            soi1=squeeze(mean(pRsq,1,'omitnan'));
            soi1=soi1(task_indices1);
            c_task_tot(1,i,channel)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
            
            part=task_indices2;
            [pRsq,~, period,coi]=wcoherence(sigPart1(1:part(end),2), sigPart2(1:part(end),2),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            soi2=squeeze(mean(pRsq,1,'omitnan'));
            soi2=soi2(task_indices2);
            c_task_tot(2,i,channel)=mean(soi2(coiidx:length(soi2)-coiidx),'omitnan');
            
            
            %~~~~~ offline complete timecourse
            
            %calc coherence
            [pRsq,~, period,coi]=wcoherence(sigPart1(:,2), sigPart2(:,2),seconds(1/fs),'PeriodLimits', [seconds(6) seconds(14)], 'VoicesPerOctave', 14);
            
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders. We'll do later.
            end
            
            %average over periods to get soi
            soi=mean(pRsq,1,'omitnan'); %soi=signal of interest, avaraged over period if interest
            
            soi_task1=soi(subtaskMat(1,1):restMat(2,1)-1);
            soi_task2=soi(subtaskMat(1,2):restMat(3,1)-1);
            soi_rest1=soi(restMat(1,1):restMat(1,2)-1);
            soi_rest2=soi(restMat(2,1):restMat(2,2)-1);
            soi_rest3=soi(restMat(3,1):restMat(3,2)-1);
            
            c_rest_tot_offmat(1, i, channel)=mean(soi_rest1(coiidx:length(soi_rest1)-coiidx),'omitnan'); %could also have a look at whole block
            c_rest_tot_offmat(2, i, channel)=mean(soi_rest2(coiidx:length(soi_rest2)-coiidx),'omitnan');
            c_rest_tot_offmat(3, i, channel)=mean(soi_rest3(coiidx:length(soi_rest3)-coiidx),'omitnan');
            c_task_tot_offmat(1,i, channel)=mean(soi_task1(coiidx:length(soi_task1)-coiidx),'omitnan');
            c_task_tot_offmat(2,i, channel)=mean(soi_task2(coiidx:length(soi_task2)-coiidx),'omitnan');
            
            
            %~~~~ count signs
            
            pos=0;
            neg=0;
            bind=0;
            
            counter=0;
            for j=1:nWinTask1
                counter=counter+1;
                diff_task_win(counter, i, channel)= c_task_win(counter,i,channel)-c_rest_tot(1,i,channel);
                if c_task_win(counter,i,channel) > c_rest_tot(1,i,channel)
                    pos=pos+1;
                elseif c_task_win(counter,i,channel)< c_rest_tot(1,i,channel)
                    neg=neg+1;
                else
                    bind=bind+1;
                end
            end
            
            for j=1:nWinTask2
                counter=counter+1;
                diff_task_win(counter, i, channel)= c_task_win(counter,i,channel)-c_rest_tot(2,i,channel);
                if c_task_win(counter,i,channel) > c_rest_tot(2,i,channel)
                    pos=pos+1;
                elseif c_task_win(counter,i,channel)< c_rest_tot(2,i,channel)
                    neg=neg+1;
                else
                    bind=bind+1;
                end
            end
            
            VZ_win(1, win, i, channel)=pos/nWinTask;
            VZ_win(2, win, i, channel)=neg/nWinTask;
            VZ_win(3, win, i, channel)=bind/nWinTask;
            
            
            %now do the same sign count for the total data online fashion
            pos=0;
            neg=0;
            bind=0;
            for j=1:2
                diff_task_tot(j,i,channel)=c_task_tot(j,i,channel)-c_rest_tot(j,i,channel);
                if c_task_tot(j,i,channel)>c_rest_tot(j,i,channel)
                    pos=pos+1;
                elseif c_task_tot(j,i,channel)<c_rest_tot(j,i,channel)
                    neg=neg+1;
                else
                    bind=bind+1;
                end
            end
            VZ_tot(1,i,channel)=pos/2;
            VZ_tot(2,i,channel)=neg/2;
            VZ_tot(3,i,channel)=bind/2;
            
            %and now the same sign count for the offline data
            pos=0;
            neg=0;
            bind=0;
            for j=1:2
                diff_task_off(j,i,channel)=c_task_tot_offmat(j,i,channel)-c_rest_tot_offmat(j,i,channel);
                if c_task_tot_offmat(j,i,channel)>c_rest_tot_offmat(j,i,channel)
                    pos=pos+1;
                elseif c_task_tot_offmat(j,i,channel)<c_rest_tot_offmat(j,i,channel)
                    neg=neg+1;
                else
                    bind=bind+1;
                end
            end
            VZ_off(1,i,channel)=pos/2;
            VZ_off(2,i,channel)=neg/2;
            VZ_off(3,i,channel)=bind/2;
            
        end %end channel
        
    end %end pairs
    
    
    %~~~~~ matlab sign test
    
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
    diffvec=reshape(diff_task_tot,[2*nPairs,nCh]);
    [p,~,stats]= signtest(diffvec(:,channel));
    results_tot(1,channel)=p; %(p,zval,sign, total);
    results_tot(2,channel)=stats.zval;
    results_tot(3,channel)=stats.sign; %number of positive signs
    results_tot(4,channel)=2*nPairs; %total number of signs
end

%total data, offline fashion
for channel=1:nCh
    diffvec=reshape(diff_task_off,[2*nPairs,nCh]);
    [p,~,stats]= signtest(diffvec(:,channel));
    results_off(1,channel)=p; %(p,zval,sign, total);
    results_off(2,channel)=stats.zval;
    results_off(3,channel)=stats.sign; %number of positive signs
    results_off(4,channel)=2*nPairs; %total number of signs
end


%% make total histogram for each channel - ratio of correct vs. incorrectly classified windows per window size

ges_results=results;
ges_results(:,7,:)=results_tot;
ges_results(:,8,:)=results_off;


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
    X=reordercats(X,windows); %apparently Matlab is reordering vector when transforming to categorical... so we need to change back.
    bar(X,ratios(:,channel),'LineWidth', 4);
    axes=gca;
    axes.FontSize=25;
    axes.FontWeight='b';
    axes.LineWidth=4;
    axes.YTick=[0.25 0.5 0.75 1];
    axes.TickLength = [0.05 0.035];
    ylim([0 1.05]);
    if chromophore==1
        title({'Ratio of correctly classified windows' , ['for HbO, channel ' num2str(channellist(channel))]}, 'FontSize', 25, 'FontWeight', 'b');
    elseif chromophore==2
        title({'Ratio of correctly classified windows', ['for HbR, channel ' num2str(channellist(channel))]}, 'FontSize', 25, 'FontWeight', 'b');
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
    writetable(mytable,['./tables/RPS_totsigntest23_channel' num2str(ch) 'HbO.csv'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/RPS_totsigntest23_channel' num2str(ch) 'HbR.csv'],'WriteRowNames',true)    
    else
        warning('Problem with chromophore identity.');
    end
  end



%% bar plot VZ - bonus plot indicating how many pairs have which ratio of correct vs. incorrect windows

fig=figure;
fig.Position(3:4)=[950,1500];

mytiles=tiledlayout(nCh*2, (length(windowsizelist)+2)/2, 'TileSpacing', 'compact');
nTiles=(length(windowsizelist)+1)*nCh;
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
        b=bar(winratio, npair, 0.5, 'LineWidth', 1.5);
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
        ylim([0 20]);
        xlabel('ratio correct feedback');
        ylabel('N dyads');
        if win==length(windowsizelist)+1
            title(['total, channel ' num2str(channellist(channel))]);
        elseif win==length(windowsizelist)+2
            title(['offline, channel ' num2str(channellist(channel))]);
        else
            title(['win ' num2str(round(windowsizelist(win)/fs)) ' s, channel ' num2str(channellist(channel))]);
        end
        
        tile=tile+1;
    end
end

if chromophore==1
    title(mytiles,'RPS: Percentage of correct classification for HbO', 'FontSize', 18, 'FontWeight', 'b');
elseif chromophore==2
    title(mytiles,'RPS: Percentage of correct classification for HbR', 'FontSize', 18, 'FontWeight', 'b');
else 
    warning('Problem with chromophore identity.');
end

linkaxes([axes{:}],'xy');