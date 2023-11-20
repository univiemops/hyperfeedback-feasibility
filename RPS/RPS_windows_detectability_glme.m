%% coherence calculation and detectability test for RPS data set.
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
% KK last version Oct. 2023


path    = ['Your/Path/To/Data']; %Unix System, change accordingly

pairlist   = {'01','02','03','04','05','06','08','09','10','11','13', ...
    '14','15','17','18','19','20','21','22','23','24','26','28','29', ...
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
channellist = [13,14,15,16];
nCh = length(channellist);

% period band of interest
poi = [6 14];

% chromophore switch
chromophore = 1; %1= HbO, 2=HbR; HbR is used in Kayhan et al. (2022)

% switch for writing out stats data and figures
writedata = true;
savefig   = false;

%initialize variables for resulting signs (VZ):
% 3: positive sign, negative + bindings
VZ_win      = NaN(3, nTestWindows, nPairs, nCh); 
VZ_tot      = NaN(3, nPairs, nCh);
VZ_off      = NaN(3, nPairs, nCh);

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
    
    load([path '/Data_' conds{co} '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test = hbo;
    load([path '/Data_' conds{co} '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test = hbo;
    
elseif chromophore==2
    load([path '/Data_' conds{co} '/RPS_01_sub1_' conds{co} '.mat']);
    hb1test = hbr;
    load([path '/Data_' conds{co} '/RPS_01_sub2_' conds{co} '.mat']);
    hb2test = hbr;
    
else
    error('chromophore not specified.');
end

%get exemplary markers - marker matrices
load([path '/Data_' conds{co} '/RPS_01_marker_' conds{co} '.mat']); 
if ~exist('t') || ~exist('fs') ||  ~exist('restMat') || ...
        ~exist('taskMat') || ~exist('subtaskMat')
    error(['failed to load in markers and/or data for pair1']);
end

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

coiidx = find(coitry > max(periodtry),1,'first'); %in samples



for win = 1:length(windowsizelist)
    
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
            
            load([path '/Data_' conds{co} '/RPS_' pairlist{pair} ...
                '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1 = hbo;
            load([path '/Data_' conds{co} '/RPS_' pairlist{pair} ...
                '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2 = hbo;
            
        elseif chromophore==2
            
            load([path '/Data_' conds{co} '/RPS_' pairlist{pair} ...
                '_sub1_' conds{co} '.mat']); %get data (hb,t)
            hb1 = hbr;
            load([path '/Data_' conds{co} '/RPS_' pairlist{pair} ...
                '_sub2_' conds{co} '.mat']); %get data (hb,t)
            hb2 = hbr;
            
        else
            error('chromophore not specified.');
        end
        
        %in case markers differ slightly for this pair, extract them again:
        load([path '/Data_' conds{co} '/RPS_' pairlist{pair} '_marker_' ...
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
            
            part = rest_indices3;
            c_rest_tot(3,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            
            %task blocks
            part = task_indices1;
            c_task_tot(1,pair,channel) = calculateCoherence(sigPart1, ...
                sigPart2, fs, poi, coiidx, part);
            
            part = task_indices2;
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
            
            %soi=signal of interest, avaraged over period if interest
            soi = mean(pRsq,1,'omitnan'); 
            
            soi_task1 = soi(subtaskMat(1,1):restMat(2,1)-1);
            soi_task2 = soi(subtaskMat(1,2):restMat(3,1)-1);
            soi_rest1 = soi(restMat(1,1):restMat(1,2)-1);
            soi_rest2 = soi(restMat(2,1):restMat(2,2)-1);
            soi_rest3 = soi(restMat(3,1):restMat(3,2)-1);
            
            c_rest_tot_offmat(1, pair, channel) = mean( ...
                soi_rest1(coiidx:length(soi_rest1)-coiidx), 'omitnan');
            c_rest_tot_offmat(2, pair, channel) = mean( ...
                soi_rest2(coiidx:length(soi_rest2)-coiidx), 'omitnan');
            c_rest_tot_offmat(3, pair, channel) = mean( ...
                soi_rest3(coiidx:length(soi_rest3)-coiidx), 'omitnan');
            c_task_tot_offmat(1, pair, channel) = mean( ...
                soi_task1(coiidx:length(soi_task1)-coiidx), 'omitnan');
            c_task_tot_offmat(2, pair, channel) = mean( ...
                soi_task2(coiidx:length(soi_task2)-coiidx), 'omitnan');
            
            
            %~~~~ count signs
            
            pos  = 0;
            neg  = 0;
            bind = 0;
            
            counter = 0;
            
            for j = 1:nWinTask1
                
                counter = counter + 1;
                diff_task_win(counter, pair, channel) = ...
                    c_task_win(counter,pair,channel) - ...
                    c_rest_tot(1,pair,channel);
                
                if c_task_win(counter,pair,channel) > ...
                        c_rest_tot(1,pair,channel)
                    pos  = pos + 1;
                    
                elseif c_task_win(counter,pair,channel) < ...
                        c_rest_tot(1,pair,channel)
                    neg  = neg + 1;
                    
                else
                    bind = bind + 1;
                end
                
            end
            
            for j = 1:nWinTask2
                
                counter = counter + 1;
                diff_task_win(counter, pair, channel) = ...
                    c_task_win(counter,pair,channel) - ...
                    c_rest_tot(2,pair,channel);
                
                if c_task_win(counter,pair,channel) > ...
                        c_rest_tot(2,pair,channel)
                    pos = pos + 1;
                    
                elseif c_task_win(counter,pair,channel) < ...
                        c_rest_tot(2,pair,channel)
                    neg = neg + 1;
                    
                else
                    bind = bind + 1;
                end
                
            end
            
            VZ_win(1, win, pair, channel) = pos/ nWinTask;
            VZ_win(2, win, pair, channel) = neg / nWinTask;
            VZ_win(3, win, pair, channel) = bind / nWinTask;
            
            
            %now do the same sign count for the total data online fashion
            pos  = 0;
            neg  = 0;
            bind = 0;
            
            for j = 1:2
                
                diff_task_tot(j,pair,channel) = ...
                    c_task_tot(j,pair,channel)-c_rest_tot(j,pair,channel);
                
                if c_task_tot(j,pair,channel) > c_rest_tot(j,pair,channel)
                    pos = pos + 1;
                    
                elseif c_task_tot(j,pair,channel) < ...
                        c_rest_tot(j,pair,channel)
                    neg = neg + 1;
                    
                else
                    bind = bind + 1;
                end
                
            end
            VZ_tot(1, pair, channel) = pos/2;
            VZ_tot(2, pair, channel) = neg/2;
            VZ_tot(3, pair, channel) = bind/2;
            
            %and now the same sign count for the offline data
            pos  = 0;
            neg  = 0;
            bind = 0;
            
            for j = 1:2
                diff_task_off(j,pair,channel) = ...
                    c_task_tot_offmat(j,pair,channel) - ...
                    c_rest_tot_offmat(j,pair,channel);
                
                if c_task_tot_offmat(j,pair,channel) > ...
                        c_rest_tot_offmat(j,pair,channel)
                    pos = pos + 1;
                    
                elseif c_task_tot_offmat(j,pair,channel) < ...
                        c_rest_tot_offmat(j,pair,channel)
                    neg = neg + 1;
                    
                else
                    bind = bind + 1;
                end
                
            end
            VZ_off(1, pair, channel) = pos/2;
            VZ_off(2, pair, channel) = neg/2;
            VZ_off(3, pair, channel) = bind/2;
            
        end %end channel
        
    end %end pairs
 
    % hierarchical binomial model; see function code below.
    teststats_glme(win,:,:) = calculateGlmeForChannels( ...
        squeeze(VZ_win(:,win,:,:)), nCh, nPairs, nWinTask);
    all_teststats_glme(win,:,:) = teststats_glme(win,:,:);

end %end window loop


% %total data, online fashion

teststats_glme_tot = calculateGlmeForChannels( ...
    VZ_tot, nCh, nPairs, 2); % nWin is two here because two runs.


% %total data, offline fashion

teststats_glme_off = calculateGlmeForChannels( ...
    VZ_off, nCh, nPairs, 2); % nWin is two here because two runs.


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
    figure;
    X=categorical(windows);
    X=reordercats(X,windows); %apparently Matlab is reordering vector when 
    %transforming to categorical... so we need to change back.
    
    bar(X,res_fig_mean(:,channel),'LineWidth', 4); hold on
    er = errorbar(X,res_fig_mean(:,channel),res_fig_std(:,channel,2), ...
        res_fig_std(:,channel,1),'LineWidth', 4);    
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
        title({'detection rate, errorbars SEM' , ['for HbO, channel ' ...
            num2str(ch)]}, 'FontSize', 25, 'FontWeight', 'b');
        
        if savefig
            saveas(gcf, ['./glme_figures/ClassAcc_HbO_' num2str(ch)], 'svg');
        else
            %do nothing
        end
    
    elseif chromophore==2
        title({'detection rate, errorbars SEM', ['for HbR, channel ' ...
            num2str(ch)]}, 'FontSize', 25, 'FontWeight', 'b');
        
        if savefig
            saveas(gcf, ['./glme_figures/ClassAcc_HbR_' num2str(ch)], 'svg');
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
    rows={'estimate'; 'SE'; 'tStat'; 'DF'; 'pValues'; 'Lower'; 'Upper'};
    mytable=array2table(squeeze(all_teststats_glme(:,channel,:))', ...
        'VariableNames', windows, 'RowNames', rows);
    
    if chromophore==1
    writetable(mytable,['./tables/RPS_glme' num2str(ch) 'HbO_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)
    elseif chromophore==2
    writetable(mytable,['./tables/RPS_glme' num2str(ch) 'HbR_' ...
        datestr(now,'mm-dd-yyyy_HH-MM') '.xlsx'],'WriteRowNames',true)    
    else
        warning('Problem with chromophore identity.');
    end
end

else 
    %do nothing
end



%% bar plot VZ - 
% supplementary plot indicating how many pairs have which ratio of correct 
% vs. incorrect windows

fig=figure;
fig.Position(3:4)=[950,1500];

mytiles=tiledlayout(nCh*2, (length(windowsizelist)+2)/2, 'TileSpacing', ...
    'compact');
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
            title(['win ' num2str(round(windowsizelist(win)/fs)) ...
                ' s, channel ' num2str(channellist(channel))]);
        end
        
        tile=tile+1;
    end
end

if chromophore==1
    title(mytiles,'RPS: Percentage of correct classification for HbO', ...
        'FontSize', 18, 'FontWeight', 'b');
elseif chromophore==2
    title(mytiles,'RPS: Percentage of correct classification for HbR', ...
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
                else
                    mat((a - 1) *nWin+b, 1) = 0;
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
