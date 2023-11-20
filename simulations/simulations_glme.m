%% simulate artificial subjects based on HRF-convolved event model plus noise
% model: rest - task - rest, rest - zero activity, task: nTrials events
% with iti intertrial interval.
% plot exemplary subjects, then simulate nIterations times nPairs of 
% subjects, calculate task windows of different sizes, determine 
% detectabililty and plot accuracy/ detection percentage.
%
% Needs SPM HRF function!
%
% KK Oct 2023 last version


%% define some basic design parameters based on Kayhan et al. 2022

nPairs=25;
poi=[6 14]; %period of interest in sec
    
nTrials=30; %number of trials
tRest=60; %60 seconds rest
fs=7.8; %sampling rate 7.8 Hz


iti=8; %8sec intertrial interval
sTask=nTrials*iti*fs; %samples task
sIti=iti*fs;
sRest=tRest*fs;
tExp=tRest+nTrials*iti+tRest;
exp_tVec=0:1/fs:tExp;

rest_sampleVec=zeros(tRest*fs,1);
sExp=length(exp_tVec);


windowsizelist=[round(50*fs) round(60*fs) round(70*fs) round(80*fs)...
    round(90*fs) round(100*fs)]; %in samples

nTestWindows=length(windowsizelist);

woffset=round(8*fs); %in samples, offset between windows

nIterations=50;


% % % % % %create basic model % % % % % %

model='peak';

modelVec=zeros(sTask+1,1);
HRF=spm_hrf(1/fs); %needs SPM


%create artifical time course via convolution
switch model
    case'peak'
        %delta peak
        for i=1:sTask
            if mod(i-1,round(sIti))==0 %15*7.8=117; mod(0,117)=0,
                %first HR at first time point
                modelVec(i)=1;
            end
        end
    case 'block'
        %block model - not investigated, use if interested 
        blockdur=4; %insert duration of each block in seconds
        boxlength=blockdur*fs;
        for i=1:sTask
            if mod(i-1,round(sIti))==0 
                for j=1:boxlength
                    if i+j>sTask+1
                        break;
                    else
                    modelVec(i+j)=1;
                    end
                end
            end
        end
end

expVec=[rest_sampleVec; modelVec; rest_sampleVec];

conv_modelVec=conv(expVec,HRF);
conv_modelVec=conv_modelVec(1:sExp); %shorten to correct length after
%convolution


% % % % % create exemplary subjects and get their basic parameters % % % %


sigPower = sum(abs(conv_modelVec(:)).^2)/numel(conv_modelVec); % linear,
%power and SNR definitions taken from Matlab function awgn
    
noiselvl=3.5;%1.55; %for scaling of gaussian noise

noisePower=(noiselvl*std(conv_modelVec))^2;
SNR=10*log10(sigPower/noisePower);
disp(SNR);


%get params from exemplary subjects
%artifical subject 1
sub1=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); 
%add random gaussian noise
%artifical subject 2
sub2=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); 
%add random gaussian noise




%%
fig=figure;
fig.Position(3:4)=[1200,1000];

tiledlayout(4,6, 'TileSpacing', 'compact');
n1=nexttile([2 4]);
plot(exp_tVec, sub1)
hold on
plot(exp_tVec, sub2)
plot(exp_tVec, conv_modelVec, 'LineWidth', 4)
xline(tRest, 'Color', 'k', 'LineWidth', 4,'LineStyle' ,'-.')
xline(tRest+nTrials*iti, 'Color', 'k','LineWidth', 4,'LineStyle' ,'-.')

hold off
xlabel('time (sec)', 'FontSize', 14, 'FontWeight', 'bold')
title('an exemplary simulated subject pair', 'FontSize', 22,...
    'FontWeight', 'b')
legend('sub1', 'sub2', 'model', 'FontSize', 15, 'Location', 'northeast')
grid on

xlim([min(exp_tVec) max(exp_tVec)]);
ylim([-0.1 0.1]);
axes=gca;
axes.FontSize=20;
axes.FontWeight='b';
axes.LineWidth=3;

%calculate basic coherence
[~,~, period,coitry]=wcoherence(sub1(1:sExp), sub2(1:sExp),...
    seconds(1/fs), 'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
    'VoicesPerOctave', 14);
%coiidx=find(coitry>max(period),1,'first'); %in samples; 
%results in same coiidx as defintion below
 
% calculate coi 
 f=1./seconds(max(period));
 cf = 6/(2*pi);
 predtimes = sqrt(2)*cf./f; % coi per wavelet frequency; in seconds
 coiidx= round(predtimes*fs); %in samples

%Values for plotting the scalogram using imagesc
[Wcoh,~, plotperiod,plotcoi]=wcoherence(sub1(1:sExp), sub2(1:sExp), ...
    seconds(1/fs), 'VoicesPerOctave', 14);

plotperiod=seconds(plotperiod); 
plotcoi=seconds(plotcoi);


%general figure:

n2=nexttile([2 2]);
    winposition=tRest+20; %in seconds
    winlength=60;
    t=exp_tVec';
    Yticks = 2.^(fix(log2(min(plotperiod))):fix(log2(max(plotperiod))));
    H=imagesc(t,log2(plotperiod),Wcoh);
    hold on
    plot(t,log2(plotcoi), 'LineWidth',5,'Color',  [0.8500 0.3250 0.0980]);
    yline(log2(6), 'LineWidth',5,'Color', 'r', 'LineStyle',':');
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.7 0.7 0.7]); 
    %fake invisible object to make rectangle appear in figure legend
    hline = line(NaN,NaN,'LineWidth',5,'Color','c'); 
    %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [winposition log2(6) winlength ...
        log2(14)-log2(6)], 'LineWidth',5,'EdgeColor','c')
    rectangle('Position', [0 log2(6) 60 log2(14)-log2(6)],...
        'LineWidth',6,'EdgeColor',[0.8 0.8 0.8])
    yline(log2(14),'LineWidth',5,'Color', 'r', 'LineStyle',':');
    xline(tRest,'LineWidth',4,'Color', 'k','LineStyle' ,'-.');
    xline(tRest+nTrials*iti,'LineWidth',4,'Color', 'k','LineStyle' ,'-.');

    hold off
    set(gca,'clim',[0 1])
    HCB=colorbar;
    set(gca,'YLim',log2([min(plotperiod),max(plotperiod)]), ...
        'YDir','reverse', 'layer','top', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Yticks'), ...
        'layer','top','FontSize',18, 'FontWeight', 'b', 'LineWidth', 3)
    ylabel('Period')
    xlabel('time (sec)');
    legend('coi for whole time course', 'period of interest',...
        'rest block', 'exemplary window', 'FontSize', 15, 'Location',...
        'northeast')
    
    title('coherence magnitude scalogram', 'FontSize', 22,...
        'FontWeight', 'b');

    
%close up figure:

n3=nexttile([1 4]);  
    t=exp_tVec';
    Yticks = 2.^(fix(log2(min(plotperiod))):fix(log2(max(plotperiod))));
    H=imagesc(t,log2(plotperiod),Wcoh);
    hold on
    %now add window and coi

    plot(t,log2(plotcoi), 'LineWidth',5,'Color', [0.8500 0.3250 0.0980]);
    
    %rest block
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.7 0.7 0.7]); 
    %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [0 log2(6) 60 log2(14)-log2(6)],...
        'LineWidth',7,'EdgeColor',[0.8 0.8 0.8])
    
    %resulting rest block
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.35 0.35 0.35]); 
    %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [coiidx/fs log2(6) 60-2*coiidx/fs...
        log2(14)-log2(6)], 'LineWidth',6,'EdgeColor', [0.35 0.35 0.35]);    
    
    %exemplary window
    hline = line(NaN,NaN,'LineWidth',5,'Color','c'); 
    %fake invisible object to make rectangle appear in figure legend    
    rectangle('Position', [winposition log2(6) winlength...
        log2(14)-log2(6)], 'LineWidth',6,'EdgeColor', 'c');

    %resulting exemplary window
    %fake invisible object to make rectangle appear in figure legend  
    hline = line(NaN,NaN,'LineWidth',5,'Color',[0.2 0.3 1]);  
    rectangle('Position', [winposition+coiidx/fs log2(6)...
        winlength-2*coiidx/fs log2(14)-log2(6)], 'LineWidth', 6,...
        'EdgeColor', [0.2 0.3 1]);
    

    plot(t+winposition,log2(plotcoi), 'LineWidth',5,'Color',...
        [0.8500 0.3250 0.0980]);   
    
    plot(t-t(end)+60,log2(plotcoi), 'LineWidth',5,'Color',...
        [0.8500 0.3250 0.0980]);
    plot(t-t(end)+winposition+winlength,log2(plotcoi),...
        'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
    xline(tRest,'LineWidth',4,'Color', 'k','LineStyle' ,'-.');
    xline(tRest+nTrials*iti,'LineWidth',4,'Color', 'k','LineStyle' ,'-.');
    yline(log2(6), 'LineWidth',8,'Color', 'r', 'LineStyle',':');
    yline(log2(14),'LineWidth',8,'Color', 'r', 'LineStyle',':');
    
    hold off
    set(gca,'clim',[0 1])
    HCB=colorbar;
  
 %    set(gca,'YLim',log2([min(plotperiod),max(plotperiod)]), ...  
    set(gca,'YLim',log2([6,14]), ...
 'YDir','reverse', 'layer','top', ...
        'YTick',log2(Yticks(:)), ...
        'YTickLabel',num2str(Yticks'), ...
        'layer','top','FontSize',20, 'FontWeight', 'b', 'LineWidth', 3)
    ylabel('Period')
    xlabel('time (sec)');
    
    title('coherence magnitude scalogram within period of interest',...
        'FontSize', 22, 'FontWeight', 'b');
legend('coi for window/block', 'rest block', 'resulting rest block',...
    'exemplary window', 'resulting window', 'FontSize', 15, ...
    'Location', 'northeast');
    
linkaxes([n1 n3], 'x')    


    
%% now simulate several measurements

%initialize variables for resulting signs (VZ):
% 3: positive sign, negative + bindings
VZ_win      = NaN(3, nTestWindows, nPairs); 
VZ_tot      = NaN(3, nPairs);
VZ_off      = NaN(3, nPairs);

%initialize for stats:
%7 outputs to glme.coefficients
teststats_glme     = NaN(nTestWindows, nIterations, 7);
all_teststats_glme = NaN(nTestWindows+2,nIterations,7);

for it=1:nIterations

% % % % % % simulate a bunch of additional subjects % % % % % 

pairMat = NaN(nPairs,2,length(sub1)); 
for p = 1:nPairs
        %sub1
        pairMat(p,1,:) = conv_modelVec+randn(size(conv_modelVec))* ...
            noiselvl*std(conv_modelVec);
        %sub2
        pairMat(p,2,:) = conv_modelVec+randn(size(conv_modelVec))* ...
            noiselvl*std(conv_modelVec);
end


%% now calc wtc windows for different window sizes

for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win); %in samples
    nWinTask=floor(sTask/(windowsize+woffset));
    
    c_task_win=NaN(nWinTask,nPairs);
    
    c_rest_tot=NaN(nPairs,1);
    c_task_tot=NaN(nPairs,1);
    
    c_rest_tot_offmat=NaN(nPairs,1);
    c_task_tot_offmat=NaN(nPairs,1);
    
    diff_task_win=NaN(nWinTask, nPairs);
    
    
    parfor p=1:nPairs
        
        % % % % % task block windows % % % % %
        
        for w=1:nWinTask
            startpoint=sRest+1+woffset; %in samples
            part=startpoint+(w-1)*woffset+(w-1)*windowsize:...
            startpoint+(w-1)*woffset+w*windowsize-1; 
            %calc coherence
            [pRsq,~, period,coi]=wcoherence(...
                squeeze(pairMat(p,1,1:part(end))),...
                squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),...
                'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
                'VoicesPerOctave', 14);
            for j=1:1:length(coi)
                pRsq(period >= coi(j), j) = NaN;
            end
            %get signal of interest
            soi=squeeze(mean(pRsq,1,'omitnan'));
            soi=soi(part);
            c_task_win(w,p)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');
        end
        
        
        
        % % % % % blocks up to end - as if online % % % % %
        
        % first rest block:
        
        part=1:sRest; %complete rest block
        [pRsq,~, period,coi]=wcoherence( ...
            squeeze(pairMat(p,1,1:part(end))), ...
            squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),...
            'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
            'VoicesPerOctave', 14);
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN;
        end
        soi1=squeeze(mean(pRsq,1,'omitnan'));
        soi1=soi1(part);
        c_rest_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        % task block:
        
        part=sRest+1:sTask;
        [pRsq,~, period,coi]=wcoherence(...
            squeeze(pairMat(p,1,1:part(end))),...
            squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),...
            'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
            'VoicesPerOctave', 14);
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN;
        end
        soi1=squeeze(mean(pRsq,1,'omitnan'));
        soi1=soi1(part);
        c_task_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        
        % % % % % for comparison: complete timecourse transformed - as
        % in offline analysis % % % % %

        [pRsq,~, period,coi] = wcoherence(squeeze(pairMat(p,1,:)),...
            squeeze(pairMat(p,2,:)),seconds(1/fs),'PeriodLimits',...
            [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
        
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN; %outer bracket
        end
        %average over periods to get soi
        soi = mean(pRsq,1,'omitnan');
        soi_rest = soi(1:sRest);
        soi_task = soi(sRest+1:sTask);
        c_rest_tot_offmat(p) = mean(soi_rest(...
            coiidx:length(soi_rest)-coiidx),'omitnan');
        c_task_tot_offmat(p) = mean(soi_task(...
            coiidx:length(soi_task)-coiidx),'omitnan');
        
    end  %end parfor pair loop        
    
    for p=1:nPairs
        
        % ~~~~ count signs

                pos=0;
                neg=0;
                bind=0;
                
                for j=1:nWinTask
                    
                    if c_task_win(j,p) > ...
                            c_rest_tot(p)
                        pos=pos+1;
                    elseif c_task_win(j,p) < ...
                            c_rest_tot(p)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                    
                end        
                
                VZ_win(1, win, p)=pos/nWinTask;
                VZ_win(2, win, p)=neg/nWinTask;
                VZ_win(3, win, p)=bind/nWinTask;
                
                
               %now do the same for the total data online fashion
               
                pos=0;
                neg=0;
                bind=0;

                    if c_task_tot(p) > c_rest_tot(p)
                        pos=pos+1;
                    elseif c_task_tot(p) < ...
                            c_rest_tot(p)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_tot(1, p)=pos;
                VZ_tot(2, p)=neg;
                VZ_tot(3, p)=bind;
                
            
              %and now the same for the offline data  
              
                pos=0;
                neg=0;
                bind=0;
                
                    if c_task_tot_offmat(p) > ...
                            c_rest_tot_offmat(p)
                        pos=pos+1;
                    elseif c_task_tot_offmat(p) < ...
                            c_rest_tot_offmat(p)
                        neg=neg+1;
                    else
                        bind=bind+1;
                    end
                VZ_off(1, p)=pos;
                VZ_off(2, p)=neg;
                VZ_off(3, p)=bind;          
        
        
    end  %end pair loop for counting signs           
        
        

    
    % hierarchical binomial model; see function code below.
    teststats_glme(win, it, :) = calculateGlme( ...
        squeeze(VZ_win(:,win,:)), nPairs, nWinTask);
    all_teststats_glme(win, it, :) = teststats_glme(win, it ,:);
    
end %end window loop

% %total data, online fashion

teststats_glme_tot = calculateGlme(VZ_tot, nPairs, 1); 


% %total data, offline fashion

teststats_glme_off = calculateGlme(VZ_off, nPairs, 1); 

all_teststats_glme(nTestWindows+1,it,:) = teststats_glme_tot;
all_teststats_glme(nTestWindows+2,it,:) = teststats_glme_off;



end %end iterations loop


%% % % % % % make histogram plot % % % % %

%initialize for plots:
res_fig_mean        = NaN(nTestWindows+2, nIterations);
res_fig_grandmean   = NaN(nTestWindows+2,1);
res_fig_std         = NaN(nTestWindows+2, nIterations ,2);
res_fig_grandstd    = NaN(nTestWindows+2,2); 
res_p               = NaN(nTestWindows+2,1);

%transform stats to inverse-logit space, then average over iterations
for win=1:nTestWindows+2
    
    for it=1:nIterations   
        
    % heights of bars:    
    res_fig_mean(win,it) = 1/ ...
        (1+exp(-(all_teststats_glme(win,it,1))));
    
    % errorbars: need to transform total (intercept plus SE) then substract
    % because matlab bar function only takes relative values.
    absError_pos = 1/( 1+exp(-(all_teststats_glme(win,it,1) + ...
        all_teststats_glme(win,it,2))) );
    
    res_fig_std(win,it,1) = absError_pos - res_fig_mean(win,it); 
    
    absError_neg = 1/( 1+exp(-(all_teststats_glme(win,it,1) - ...
        all_teststats_glme(win,it,2))) );
    
    res_fig_std(win,it,2) = res_fig_mean(win,it) - absError_neg;     
    
    end
    
    res_p(win) = nanmean(all_teststats_glme(win,:,7));
    
    res_fig_grandmean(win)  = nanmean(res_fig_mean(win,:));
    res_fig_grandstd(win,1) = nanmean(res_fig_std(win,:,1));
    res_fig_grandstd(win,2) = nanmean(res_fig_std(win,:,2));
    
end

disp(res_p);

% all bars x axis 
xratios=zeros(size(all_teststats_glme,1),1);

%make a figure for every channel
    
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
    
    bar(X,res_fig_grandmean(:),'LineWidth', 4); hold on
    er = errorbar(X,res_fig_grandmean(:),res_fig_grandstd(:,2), ...
        res_fig_grandstd(:,1),'LineWidth', 4);    
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

    title({'classification accuracy for 25 sim. dyads',...
        '50 iterations'}, 'FontSize', 25, 'FontWeight', 'b');
    legend('average accuracy', 'average SEM','FontSize', 15,...
        'Location','southeast');





%% function definitions

% mixed model
function teststats_glme = calculateGlme(VZ_mat, nPairs, ...
    nWin)

    teststats_glme = NaN(7,1);
        
        % Extract y values from VZ_tot
        y = squeeze(VZ_mat(1, :))*nWin;

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
        teststats_glme(:) = table2array(conv(1, 2:8));
        
        %for the simulations with low SNR and accuracy rate of almost
        %always one, sometimes the solvers fail even though a 1-2 values
        %are not 1. I will remove these results for now.
         if teststats_glme(5)>0.9 
             teststats_glme(:) = NaN(7,1);
        end
end
