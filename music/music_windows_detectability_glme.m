%% coherence calculation and detectability test for the music data set.
% This script calculates the wavelet tranform coherence within a certain 
%frequency band (period of interest) for different task windows.
% It then calculates the difference in coherence between the task
% window/unit of interest and the complete rest block and determines the 
% sign (VZ=Vorzeichen) of the coherence difference. Positive sign is a 
% success, negative a fail, as would be the case in a binary feedback.
% Statistics are assessed via a generalized linear mixed model under the
% null hypothesis that the detection accuracy is at chance.
% The accuracy/ detection percentage are plotted for illustration.
% 
% last version KK Oct 2023



path=['/Your/Path/To/Data']; 
%use your data path; Unix System used, adapt accordingly

%experiment specifics
pairlist={'01','02','03','04','05','06', '07','08','09','10','11','12'};
nPairs=length(pairlist);

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

%initialize variables for results of sign test 
%(coherence_task > coherence_rest) and resulting signs (VZ)
results=NaN(4,nTestWindows,nCh); %(p,zval,sign, total);
results_tot=NaN(4,nCh);
results_off=NaN(4,nCh);
VZ_win=NaN(3, nTestWindows, nPairs, nCh); %3: positive sign, 
%negative + bindings
VZ_tot=NaN(3,nPairs,nCh);
VZ_off=NaN(3,nPairs,nCh);

%teststats_glm=NaN(nTestWindows,nCh,4); %4 outputs to mdl.coefficients: estimate, SE, tstat, pValue
teststats_glme=NaN(nTestWindows,nCh,7); %7 outputs to glme.coefficients
all_teststats_glme = NaN(nTestWindows+2, nCh, 7);

%initialize for plots:
res_fig_mean = NaN(nTestWindows+2, nCh);
res_fig_std  = NaN(nTestWindows+2, nCh,2);

%nWin_task_vec=NaN(nTestWindows,1);

%chromophore switch
chromophore = 1; %1= HbO, 2=HbR
%in paper HbO is used

writedata = true;
savefig   = false;
    
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

coiidx=find(coitry>max(periodtry),1,'first'); %in samples

%%
for win=1:length(windowsizelist)
    
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
    diff_task_tot=NaN(nPairs,nCh);
    diff_task_off=NaN(nPairs,nCh);
    
    
    
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
               startpoint = 1800+200+woffset; %in samples
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

                %average over periods to get soi
                soi = mean(pRsq,1,'omitnan'); %soi=signal of interest

                soi_task1 = soi(task_indices);
                soi_rest1 = soi(rest_indices);

                c_rest_tot_offmat(pair, channel ) = mean( ...
                    soi_rest1(coiidx:length(soi_rest1)-coiidx),'omitnan'); 
                c_task_tot_offmat(pair, channel) = mean( ...
                    soi_task1(coiidx:length(soi_task1)-coiidx),'omitnan');

            
       
        % ~~~~ count signs

                pos=0;
                neg=0;
                bind=0;
                
                counter=0;
                
                for j=1:nWinTask
                    
                    counter=counter+1;
                    
                    diff_task_win(counter, pair, channel) = ...
                        c_task_win(j,pair,channel) - ...
                        c_rest_tot(pair,channel);
                    if c_task_win(j,pair,channel) > ...
                            c_rest_tot(pair,channel)
                        pos=pos+1;
                    elseif c_task_win(j,pair,channel) < ...
                            c_rest_tot(pair,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                    
                end        
                
                VZ_win(1, win, pair, channel)=pos/nWinTask;
                VZ_win(2, win, pair, channel)=neg/nWinTask;
                VZ_win(3, win, pair, channel)=bind/nWinTask;
                
                
                
               %now do the same for the total data online fashion
               
                pos=0;
                neg=0;
                bind=0;
               
                    diff_task_tot(pair,channel) = ...
                        c_task_tot(pair,channel) - ...
                        c_rest_tot(pair,channel);
                    
                    if c_task_tot(pair,channel) > c_rest_tot(pair,channel)
                        pos=pos+1;
                    elseif c_task_tot(pair,channel) < ...
                            c_rest_tot(pair,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_tot(1, pair, channel)=pos;
                VZ_tot(2, pair, channel)=neg;
                VZ_tot(3, pair, channel)=bind;
                
            
                
              %and now the same for the offline data  
              
                pos=0;
                neg=0;
                bind=0;
                
                    diff_task_off(pair,channel) = ...
                        c_task_tot_offmat(pair,channel) - ...
                        c_rest_tot_offmat(pair,channel);
                    if c_task_tot_offmat(pair,channel) > ...
                            c_rest_tot_offmat(pair,channel)
                        pos=pos+1;
                    elseif c_task_tot_offmat(pair,channel) < ...
                            c_rest_tot_offmat(pair,channel)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_off(1, pair, channel)=pos;
                VZ_off(2, pair, channel)=neg;
                VZ_off(3, pair, channel)=bind;  
                
                
                
        
        end %end channelloop
        
    end%end pair loop
     
    % hierarchical binomial model; see function code below.
    teststats_glme(win,:,:) = calculateGlmeForChannels( ...
        squeeze(VZ_win(:,win,:,:)), nCh, nPairs, nWinTask);
    all_teststats_glme(win,:,:) = teststats_glme(win,:,:);

end %end window loop

% %total data, online fashion

teststats_glme_tot = calculateGlmeForChannels( ...
    VZ_tot, nCh, nPairs, 1); % nWin is one


% %total data, offline fashion

teststats_glme_off = calculateGlmeForChannels( ...
    VZ_off, nCh, nPairs, 1); % nWin is one

all_teststats_glme(nTestWindows+1,:,:) = teststats_glme_tot;
all_teststats_glme(nTestWindows+2,:,:) = teststats_glme_off;

  
%% plot with error bars
% make total histogram for each channel - ratio of correct vs. incorrectly
% classified windows per window size

% for bar plots get inverse logit links of glme results.

for win=1:nTestWindows+2
    for channel=1:nCh   
        
    % heights of bars:    
    res_fig_mean(win,channel) = 1/ ...
        (1+exp(-(all_teststats_glme(win,channel,1))));
    
    % errorbars: need to transform total (intercept plus SE) then substract
    % because matlab bar function only takes relative values.
    absError_pos = 1/( 1+exp(-(all_teststats_glme(win,channel,1) + ...
        all_teststats_glme(win,channel,2))) );
    res_fig_std(win,channel,1) = absError_pos - res_fig_mean(win,channel);   
    absError_neg = 1/( 1+exp(-(all_teststats_glme(win,channel,1) - ...
        all_teststats_glme(win,channel,2))) );
    res_fig_std(win,channel,2) = res_fig_mean(win,channel) - absError_neg;     
    
    end
end


% all bars x axis 
xratios=zeros(size(all_teststats_glme,1),nCh);

%make a figure for every channel
for ch=1:nCh
    
    for i=1:size(all_teststats_glme,1)

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
    X=reordercats(X,windows); %apparently Matlab is reordering vector when 
    %transforming to categorical... so we need to change back.
    
    bar(X,res_fig_mean(:,ch),'LineWidth', 4); hold on
    er = errorbar(X,res_fig_mean(:,ch),res_fig_std(:,ch,2), ...
        res_fig_std(:,ch,1),'LineWidth', 4);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    yline(0,'LineWidth', 4); 
    hold off
    
    axes=gca;
    axes.FontSize=25;
    axes.FontWeight='b';
    axes.LineWidth=4;
    axes.YTick=[0.25 0.5 0.75 1];
    yline(0.5,'LineWidth', 2);
  %  axes.TickLength = [0.05 0.035];
    ylim([0 1.05]);
    
    if chromophore==1
        title({'classification accuracy, errorbars SEM' , ['for HbO, channel ' ...
            num2str(channellist(ch))]}, 'FontSize', 25, 'FontWeight', 'b');
        
        if savefig
            saveas(gcf, ['./figures/ClassAcc_HbO_' num2str(ch)], 'svg');
        else
            %do nothing
        end
        
    elseif chromophore==2
        title({'classification accuracy, errorbars SEM', ['for HbR, channel ' ...
            num2str(channellist(ch))]}, 'FontSize', 25, 'FontWeight', 'b');
        
        if savefig
            saveas(gcf, ['./figures/ClassAcc_HbR_' num2str(ch)], 'svg');
        else
            %do nothing
        end
        
    else
        warning('Problem with chromophore identity.');
    end
    
    
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
    rows={'estimate'; 'SE'; 'tStat'; 'DF'; 'pValues'; 'Lower'; 'Upper'};
    mytable=array2table(squeeze(all_teststats_glme(:,channel,:))', ...
        'VariableNames', windows, 'RowNames', rows);
    
    if chromophore==1
    writetable(mytable,['./tables/music_glme_contcorr_' num2str(ch) 'HbO_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/music_glme_contcorr_' num2str(ch) 'HbR_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)   
    else
        warning('Problem with chromophore identity.');
    end
end

else
    %do nothing
end


%% bar plot VZ - bonus plot indicating how many pairs have which ratio of 
%correct vs. incorrect data

fig=figure;
%fig.Position(3:4)=[8000,600];
fig.Position(3:4)=[950,1500];

%mytiles=tiledlayout(nCh, length(windowsizelist)+2, 'TileSpacing', 'compact');
mytiles=tiledlayout(nCh*2, (length(windowsizelist)+2)/2, 'TileSpacing', ...
    'compact');
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
   title(['win ' num2str(round(windowsizelist(win)/fs)) ' s, channel ' ...
       num2str(channellist(channel))]);
end

xlabel('ratio correct feedback');
ylabel('N dyads');
tile=tile+1;
end
end

if chromophore==1
    title(mytiles,'Music: Percentage of correct classification for HbO',...
        'FontSize', 18, 'FontWeight', 'b');
elseif chromophore==2
    title(mytiles,'Music: Percentage of correct classification for HbR',...
        'FontSize', 18, 'FontWeight', 'b');
else 
    warning('Problem with chromophore identity.');
end

linkaxes([axes{:}],'xy');

%% function definitions

% coherence
function c_value = calculateCoherence(sigPart1, sigPart2, fs, poi, ...
    coiidx, part)

    % Calculate coherence with Matlab's coherence function
    [pRsq, ~, period, coi] = wcoherence(sigPart1(1:part(end), 2), ...
        sigPart2(1:part(end), 2), seconds(1/fs), 'PeriodLimits', ...
        [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);

    % NaN everything outside of coi as safety net
    for j = 1:1:length(coi)
        pRsq(period >= coi(j), j) = NaN;
    end

    % Get signal of interest = soi
    soi = squeeze(mean(pRsq, 1, 'omitnan'));
    soi = soi(part);

    c_value = mean(soi(coiidx:length(soi) - coiidx), 'omitnan');
end


% mixed model
function teststats_glme = calculateGlmeForChannels(VZ_mat, nCh, nPairs, ...
    nWin)

    teststats_glme = NaN(nCh, 7);

    for channel = 1:nCh
        
        % Extract y values from VZ_tot
        y = squeeze(VZ_mat(1, :, channel))*nWin;

        % Initialize the result matrix for this channel
        mat = NaN(nPairs * nWin, 2);

        % Populate the result matrix
        for a = 1:nPairs
            for b = 1:nWin
                
                if y(a) >= b
                    mat((a - 1) *nWin + b, 1) = 1;
                elseif y(a) < b
                    mat((a - 1) *nWin + b, 1) = 0;
                else %if NaN etc
                    mat((a - 1) *nWin + b, 1) = NaN;
                end
                
                mat((a - 1) *nWin+b, 2) = a;
            end
        end
        
        %continuity correction in case all windows are successes or
        %failures: add two elements to last pair (doesn't matter which one)
        %one success and one failure
        if ~any(mat(:,1)>0) || ~any(mat(:,1)<1) 
            addme = [1, nPairs; 0, nPairs];
            mat = [mat; addme];
        else
            %do nothing, business as usual
        end

        % Create a table from the result matrix
        tab = array2table(mat, 'VariableNames', {'success', 'dyadID'});

        % Fit a generalized linear mixed model
        glme = fitglme(tab, 'success ~ 1 + (1|dyadID)', 'Distribution', ...
            'Binomial');
        
        % Display the model
        % disp(glme);

        % Due to data set format, we have to convert first to table, then
        % array
        conv = dataset2table(glme.Coefficients);
        teststats_glme(channel, :) = table2array(conv(1, 2:8));
        
    end
end
