%% simulate artificial subjects based on HRF-convolved event model plus noise
% model: rest - task - rest, rest - zero activity, task: nTrials events
% with iti intertrial interval.
% plot exemplary subjects, then simulate nIterations times nPairs of 
% subjects, calculate task windows of different sizes, calculate mixed
% model, make bar plot.
% Needs SPM HRF function!
%
%
% KK 2024


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


woffset=round(8*fs); %in samples, offset between windows

nTestWindows=length(windowsizelist);
nWinTaskMax=floor(sTask/(windowsizelist(1)+woffset));

nIterations=50;

datacell{1}=NaN(nWinTaskMax*nPairs,nTestWindows+2);
datacell{2}=NaN(nWinTaskMax*nPairs,nTestWindows+2);

bigmat = NaN(nWinTaskMax*nPairs, 3, nTestWindows+2);

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
    
noiselvl=1.55;%6.25;%5;%1.55;3.5; %for scaling of gaussian noise

noisePower=(noiselvl*std(conv_modelVec))^2;
SNR=10*log10(sigPower/noisePower);
disp(SNR);

%%
%get params from exemplary subjects
%artifical subject 1
sub1=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); 
%add random gaussian noise
%artifical subject 2
sub2=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); 
%add random gaussian noise



%calculate basic coherence
[~,~, period,coitry]=wcoherence(sub1(1:sExp), sub2(1:sExp),...
    seconds(1/fs), 'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
    'VoicesPerOctave', 14);
%coiidx=find(coitry>max(period),1,'first'); %in samples; 
%results is the same coiidx as defintion below
 
% calculate coi 
 f=1./seconds(max(period));
 cf = 6/(2*pi);
 predtimes = sqrt(2)*cf./f; % coi per wavelet frequency; in seconds
 coiidx= round(predtimes*fs); %in samples


    
%% now simulate several measurements

%initialize for stats:
%7 outputs to glme.coefficients
%for continuous stats
teststats_glme_cont     = NaN(nTestWindows, nIterations, 7);
all_teststats_glme_cont = NaN(nTestWindows+2,nIterations,7);

%for classification stats
teststats_glme     = NaN(nTestWindows, nIterations, 7);
all_teststats_glme = NaN(nTestWindows+2,nIterations,7);

    % % % % % % simulate a bunch of additional subjects % % % % %
    
   s=rng; %set seed
 %  rng(s);
   %
    pairMat = NaN(nPairs,2,length(sub1),nIterations);
    
for l=1:nIterations
    

        
    for p = 1:nPairs
        %sub1
        pairMat(p,1,:, l) =conv_modelVec +randn(size(conv_modelVec))* ...
            noiselvl*std(conv_modelVec);
        %sub2
        pairMat(p,2,:, l) =conv_modelVec + randn(size(conv_modelVec))* ...
            noiselvl*std(conv_modelVec);
    end
    
end


%% now calc wtc windows for different window sizes
for it=1:nIterations
    
for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win); %in samples
    nWinTask=floor(sTask/(windowsize+woffset));
    
    c_task_win=NaN(nWinTask,nPairs);
    
    c_rest_tot=NaN(nPairs,1);
    c_task_tot=NaN(nPairs,1);
    
    c_rest_tot_offmat=NaN(nPairs,1);
    c_task_tot_offmat=NaN(nPairs,1);
    
    diff_task_win=NaN(nWinTask, nPairs);
    diff_task_tot = NaN(1,nPairs);
    diff_task_off = NaN(1,nPairs);
    
    VZ_task_win      = NaN(nWinTask, nPairs); 
    VZ_task_tot      = NaN(1,nPairs);
    VZ_task_off      = NaN(1,nPairs);
    
    for p=1:nPairs
        
        % % % % % task block windows % % % % %
        
        for w=1:nWinTask
            startpoint=sRest+1+woffset; %in samples
            part=startpoint+(w-1)*woffset+(w-1)*windowsize:...
            startpoint+(w-1)*woffset+w*windowsize-1; 
            %calc coherence
            [pRsq,~, period,coi]=wcoherence(...
                squeeze(pairMat(p,1,part, it)),...
                squeeze(pairMat(p,2,part, it)),seconds(1/fs),...
                'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
                'VoicesPerOctave', 14);
%             for j=1:1:length(coi)
%                 pRsq(period >= coi(j), j) = NaN;
%             end
                        
            %reduced pRsq to avoid mayer wave freq
            redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
            
            
            %get signal of interest
            soi1=squeeze(mean(redpRsq,1,'omitnan'));

             
            soi=soi1;%(part);
            c_task_win(w,p)=mean(soi(coiidx:length(soi)-coiidx),'omitnan');

         title(['windowsize: ' num2str(round(windowsize/fs)) ' mean: ' num2str(c_task_win(w,p))]);
        end
        
        
        
        % % % % % blocks up to end - as if online % % % % %
        
        % first rest block:
        
        part=1:sRest; %complete rest block
        [pRsq,~, period,coi]=wcoherence( ...
            squeeze(pairMat(p,1,part, it)), ...
            squeeze(pairMat(p,2,part, it)),seconds(1/fs),...
            'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
             'VoicesPerOctave', 14);
%         for j=1:1:length(coi)
%             pRsq(period >= coi(j), j) = NaN;
%         end
        
        %reduced pRsq to avoid mayer wave freq
        redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
    
        soi1=squeeze(mean(redpRsq,1,'omitnan'));
        soi1=soi1;%(part);
        c_rest_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        % task block:
        
        part=sRest+1:sRest+sTask;
        [pRsq,~, period,coi]=wcoherence(...
            squeeze(pairMat(p,1,part, it)),...
            squeeze(pairMat(p,2,part, it)),seconds(1/fs),...
            'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],...
            'VoicesPerOctave', 14);
%         for j=1:1:length(coi)
%             pRsq(period >= coi(j), j) = NaN;
%         end
        
        %reduced pRsq to avoid mayer wave freq
        redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
        
        soi1=squeeze(mean(redpRsq,1,'omitnan'));
        soi1=soi1;%(part);
        c_task_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        
        % % % % % for comparison: complete timecourse transformed - as
        % in offline analysis % % % % %

        [pRsq,~, period,coi] = wcoherence(squeeze(pairMat(p,1,:, it)),...
            squeeze(pairMat(p,2,:, it)),seconds(1/fs),'PeriodLimits',...
            [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
        
%         for j=1:1:length(coi)
%             pRsq(period >= coi(j), j) = NaN; %outer bracket
%         end

        %reduced pRsq to avoid mayer wave freq
        redpRsq = pRsq((seconds(period)<=9 | seconds(period)>=11), :);
        
        %average over periods to get soi
        soi = mean(redpRsq,1,'omitnan');
%         soi_rest = soi(1:sRest);
%         soi_task = soi(sRest+1:sTask);
        soi_rest = soi(1:sRest);
        soi_task = soi(sRest+1:sRest+sTask);
        c_rest_tot_offmat(p) = mean(soi_rest(...
            coiidx:length(soi_rest)-coiidx),'omitnan');
        c_task_tot_offmat(p) = mean(soi_task(...
            coiidx:length(soi_task)-coiidx),'omitnan');
        
    end  %end parfor pair loop        
    
    for p=1:nPairs
        
                
                for j=1:nWinTask
                    
                    if c_task_win(j,p) > ...
                            c_rest_tot(p)
                    elseif c_task_win(j,p) < ...
                            c_rest_tot(p)
                    else
                        %do nothing
                    end
                 diff_task_win(j,p) = c_task_win(j,p) - c_rest_tot(p);                   
                end        
                


                 %now do the same for the total data online fashion

                    if c_task_tot(p) > c_rest_tot(p)
                    elseif c_task_tot(p) < ...
                            c_rest_tot(p)                     
                    else
                        %do nothing
                    end
                
                diff_task_tot(p) = c_task_tot(p) - c_rest_tot(p);                
            
              %and now the same for the offline data  

                
                    if c_task_tot_offmat(p) > ...
                            c_rest_tot_offmat(p)
 
                    elseif c_task_tot_offmat(p) < ...
                            c_rest_tot_offmat(p)
                    else
                        %do nothing
                    end
                
                diff_task_off(p) = c_task_tot_offmat(p) - c_rest_tot_offmat(p);
        
    end  %end pair loop          
 

    [teststats_glme_cont(win, it, :), contWinTab] = calculateGlmeCont( ...
        diff_task_win, nPairs, nWinTask);
    all_teststats_glme_cont(win, it, :) = teststats_glme_cont(win, it ,:);
    

    
    
end %end window loop

% % %total data, online fashion
[teststats_glme_tot_cont, contWinTabtot] = calculateGlmeCont(diff_task_tot, nPairs, 1); 
    

% %total data, offline fashion 
[teststats_glme_off_cont, contWinTaboff] = calculateGlmeCont(diff_task_off, nPairs, 1); 
    

all_teststats_glme_cont(nTestWindows+1,it,:) = teststats_glme_tot_cont;
all_teststats_glme_cont(nTestWindows+2,it,:) = teststats_glme_off_cont;





end %end iterations loop

%%
%save('workspace_SNR7.mat')



%% % % % % % make histogram plot - continuous % % % % %

%initialize for plots:
res_fig_mean_cont        = NaN(nTestWindows+2, nIterations);
res_fig_grandmean_cont   = NaN(nTestWindows+2,1);
res_fig_std_cont         = NaN(nTestWindows+2, nIterations ,2);
res_fig_grandstd_cont    = NaN(nTestWindows+2,2); 
res_p_cont               = NaN(nTestWindows+2,1);


for win=1:nTestWindows+2
    

    for it=1:nIterations  
        
    % heights of bars:    
    res_fig_mean_cont(win,it) = (all_teststats_glme_cont(win,it,1));

    res_fig_std_cont(win,it,1) = all_teststats_glme_cont(win,it,2);
    res_fig_std_cont(win,it,2) = all_teststats_glme_cont(win,it,2);    
    
    end
    
    res_p_cont(win) = nanmean(all_teststats_glme_cont(win,:,5));
    
    res_fig_grandmean_cont(win)  = nanmean(res_fig_mean_cont(win,:));
    res_fig_grandstd_cont(win,1) = nanmean(res_fig_std_cont(win,:,1));
    res_fig_grandstd_cont(win,2) = nanmean(res_fig_std_cont(win,:,2));
    
end

disp(res_p_cont);

% all bars x axis 
xratios=zeros(size(all_teststats_glme_cont,1),1);

    
    for i=1:size(all_teststats_glme_cont,1)

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
    
    bar(X,res_fig_grandmean_cont(:),'LineWidth', 4); hold on
    er = errorbar(X,res_fig_grandmean_cont(:),res_fig_grandstd_cont(:,2), ...
        res_fig_grandstd_cont(:,1),'LineWidth', 4);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';  
    yline(0,'LineWidth', 4); 
    hold off
    
    axes=gca;
    axes.FontSize=25;
    axes.FontWeight='b';
    axes.LineWidth=4;
   % axes.YTick=[0.1 0.2];
    axes.YTick=[0.1 0.2 0.3 0.4 0.5 0.6];
%    yline(0.5,'LineWidth', 2);
  %  axes.TickLength = [0.05 0.035];
    ylim([0 1]);

    title({'\Delta WTC for 25 dyads',['SNR=' num2str(round(SNR)) ', ' num2str(nIterations) ' sim. data sets']},'FontSize', 25, 'FontWeight', 'b');
    legend(' \Delta WTC', 'SEM','FontSize', 15,...
        'Location','southeast');
    
 writedata=true;
    
 
 
 if writedata
     
    data = [res_fig_grandmean_cont, squeeze(res_fig_grandstd_cont(:,1)), res_p_cont];
        rows={'\Delta WTC'; 'SE';'pValues'};
    mytable=array2table(data', ...
        'VariableNames', windows, 'RowNames', rows); 
        writetable(mytable,['./tables/sim_' num2str(round(SNR,2)) '_' datestr(now,'mm-dd-yyyy_HH-MM') '.csv'],'WriteRowNames',true) 
 end




%% function definitions



% mixed model - continuous
function [teststats_glme_cont, outtab] = calculateGlmeCont(diff_task_win, nPairs, ...
    nWin)

    teststats_glme_cont = NaN(7,1);
        
        % Extract y values from VZ_tot
        y = diff_task_win;

        % Initialize the result matrix for this channel
        mat = NaN(nPairs * nWin, 2);

        % Populate the result matrix
       for a = 1:nPairs
            for b = 1:nWin
                
                mat((a - 1) * nWin + b , 1) = atanh(y(b, a));                

                mat((a - 1) *nWin+b, 2) = a;
            end
        end
        % Create a table from the result matrix
        tab = array2table(mat, 'VariableNames', {'DeltaWTC', 'dyadID'});
        outtab = tab(:,[2,1]);

        % Fit a generalized linear mixed model
        glme = fitglme(tab, 'DeltaWTC ~ 1 + (1|dyadID)', 'Distribution', ...
            'normal');
        
        % Display the model
        % disp(glme);

        % Due to data set format, we have to convert first to table, then
        % array
        conv = dataset2table(glme.Coefficients);
        teststats_glme_cont(:) = table2array(conv(1, 2:8));
        

end

