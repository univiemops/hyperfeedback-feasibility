%% simulate artificial subjects based on HRF-convolved event model plus noise
% model: rest - task - rest, rest - zero activity, task: nTrials events
% with iti intertrial interval.
% plot exemplary subjects, simulate nIterations times nPairs of subjects,
% calculate task windows of different sizes, plot histogram of ratio of
% correct detection/feedback.
% Needs SPM HRF function!
%
% KK May 2023


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


windowsizelist=[round(50*fs) round(60*fs) round(70*fs) round(80*fs) round(90*fs) round(100*fs)]; %in samples

nTestWindows=length(windowsizelist);

woffset=round(8*fs); %in samples, offset between windows

nIterations=50;

results=NaN(4,nTestWindows,nIterations); %(p,zval,sign, total);
results_tot=NaN(4,nIterations);
results_off=NaN(4,nIterations);
VZ_win=NaN(3, nTestWindows, nPairs, nIterations); %3: positive sign, negative + bindings
VZ_tot=NaN(3,nPairs, nIterations);
VZ_off=NaN(3,nPairs, nIterations);



% % % % % %create basic model % % % % % %

model='peak';

modelVec=zeros(sTask+1,1);
HRF=spm_hrf(1/fs); %needs SPM


%create artifical time course via convolution
switch model
    case'peak'
        %delta peak
        for i=1:sTask
            if mod(i-1,round(sIti))==0 %15*7.8=117; mod(0,117)=0, first HR at first time point
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
conv_modelVec=conv_modelVec(1:sExp); %shorten to correct length after convolution


% % % % % create exemplary subjects and get their basic parameters % % % %


sigPower = sum(abs(conv_modelVec(:)).^2)/numel(conv_modelVec); % linear, power and SNR definitions taken from Matlab function awgn
    
noiselvl=3.5; %for scaling of gaussian noise

noisePower=(noiselvl*std(conv_modelVec))^2;
SNR=10*log10(sigPower/noisePower);
disp(SNR);


%get params from exemplary subjects
%artifical subject 1
sub1=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); %add random gaussian noise
%artifical subject 2
sub2=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec); %add random gaussian noise




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
title('an exemplary simulated subject pair', 'FontSize', 22, 'FontWeight', 'b')
legend('sub1', 'sub2', 'model', 'FontSize', 15, 'Location', 'northeast')
grid on

xlim([min(exp_tVec) max(exp_tVec)]);
ylim([-0.1 0.1]);
axes=gca;
axes.FontSize=20;
axes.FontWeight='b';
axes.LineWidth=3;

%calculate basic coherence
[~,~, period,coitry]=wcoherence(sub1(1:sExp), sub2(1:sExp),seconds(1/fs), 'PeriodLimits', [seconds(poi(1)) seconds(poi(2))], 'VoicesPerOctave', 14);
%coiidx=find(coitry>max(period),1,'first'); %in samples; results in same coiidx as defintion below
 
% calculate coi 
 f=1./seconds(max(period));
 cf = 6/(2*pi);
 predtimes = sqrt(2)*cf./f; % coi per wavelet frequency; in seconds
 coiidx= round(predtimes*fs); %in samples

%Values for plotting the scalogram using imagesc
[Wcoh,~, plotperiod,plotcoi]=wcoherence(sub1(1:sExp), sub2(1:sExp),seconds(1/fs), 'VoicesPerOctave', 14);

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
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.7 0.7 0.7]); %fake invisible object to make rectangle appear in figure legend
    hline = line(NaN,NaN,'LineWidth',5,'Color','c'); %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [winposition log2(6) winlength log2(14)-log2(6)], 'LineWidth',5,'EdgeColor','c')
    rectangle('Position', [0 log2(6) 60 log2(14)-log2(6)], 'LineWidth',6,'EdgeColor',[0.8 0.8 0.8])
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
    legend('coi for whole time course', 'period of interest', 'rest block', 'exemplary window', 'FontSize', 15, 'Location', 'northeast')
    
    title('coherence magnitude scalogram', 'FontSize', 22, 'FontWeight', 'b');

    
%close up figure:

n3=nexttile([1 4]);  
    t=exp_tVec';
    Yticks = 2.^(fix(log2(min(plotperiod))):fix(log2(max(plotperiod))));
    H=imagesc(t,log2(plotperiod),Wcoh);
    hold on
    %now add window and coi

    plot(t,log2(plotcoi), 'LineWidth',5,'Color', [0.8500 0.3250 0.0980]);
    
    %rest block
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.7 0.7 0.7]); %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [0 log2(6) 60 log2(14)-log2(6)], 'LineWidth',7,'EdgeColor',[0.8 0.8 0.8])
    
    %resulting rest block
    hline = line(NaN,NaN,'LineWidth',5,'Color',	[0.35 0.35 0.35]); %fake invisible object to make rectangle appear in figure legend
    rectangle('Position', [coiidx/fs log2(6) 60-2*coiidx/fs log2(14)-log2(6)], 'LineWidth',6,'EdgeColor', [0.35 0.35 0.35]);    
    
    %exemplary window
    hline = line(NaN,NaN,'LineWidth',5,'Color','c'); %fake invisible object to make rectangle appear in figure legend    
    rectangle('Position', [winposition log2(6) winlength log2(14)-log2(6)], 'LineWidth',6,'EdgeColor', 'c');

    %resulting exemplary window
    hline = line(NaN,NaN,'LineWidth',5,'Color',[0.2 0.3 1]); %fake invisible object to make rectangle appear in figure legend   
    rectangle('Position', [winposition+coiidx/fs log2(6) winlength-2*coiidx/fs log2(14)-log2(6)], 'LineWidth',6,'EdgeColor', [0.2 0.3 1]);
    

    plot(t+winposition,log2(plotcoi), 'LineWidth',5,'Color', [0.8500 0.3250 0.0980]);   
    
    plot(t-t(end)+60,log2(plotcoi), 'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
    plot(t-t(end)+winposition+winlength,log2(plotcoi), 'LineWidth',5,'Color',[0.8500 0.3250 0.0980]);
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
    
    title('coherence magnitude scalogram within period of interest', 'FontSize', 22, 'FontWeight', 'b');
legend('coi for window/block', 'rest block', 'resulting rest block', 'exemplary window', 'resulting window', 'FontSize', 15,'Location', 'northeast');
    
linkaxes([n1 n3], 'x')    


    
%% now simulate several measurements

for it=1:nIterations

% % % % % % simulate a bunch of additional subjects % % % % % 

pairMat=NaN(nPairs,2,length(sub1)); 
for p=1:nPairs
        %sub1
        pairMat(p,1,:)=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec);
        %sub2
        pairMat(p,2,:)=conv_modelVec+randn(size(conv_modelVec))*noiselvl*std(conv_modelVec);
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
            part=startpoint+(w-1)*woffset+(w-1)*windowsize:startpoint+(w-1)*woffset+w*windowsize-1; 
            %calc coherence
            [pRsq,~, period,coi]=wcoherence(squeeze(pairMat(p,1,1:part(end))), squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
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
        [pRsq,~, period,coi]=wcoherence(squeeze(pairMat(p,1,1:part(end))), squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN;
        end
        soi1=squeeze(mean(pRsq,1,'omitnan'));
        soi1=soi1(part);
        c_rest_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        % task block:
        
        part=sRest+1:sTask;
        [pRsq,~, period,coi]=wcoherence(squeeze(pairMat(p,1,1:part(end))), squeeze(pairMat(p,2,1:part(end))),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN;
        end
        soi1=squeeze(mean(pRsq,1,'omitnan'));
        soi1=soi1(part);
        c_task_tot(p)=mean(soi1(coiidx:length(soi1)-coiidx),'omitnan');
        
        
        % % % % % for comparison: complete timecourse transformed - as in offline analysis % % % % %
        
        [pRsq,~, period,coi]=wcoherence(squeeze(pairMat(p,1,:)), squeeze(pairMat(p,2,:)),seconds(1/fs),'PeriodLimits', [seconds(poi(1)) seconds(poi(2))],'VoicesPerOctave', 14);
        
        for j=1:1:length(coi)
            pRsq(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders. We'll do later.
        end
        %average over periods to get soi
        soi=mean(pRsq,1,'omitnan');
        soi_rest=soi(1:sRest);
        soi_task=soi(sRest+1:sTask);
        c_rest_tot_offmat(p)=mean(soi_rest(coiidx:length(soi_rest)-coiidx),'omitnan');
        c_task_tot_offmat(p)=mean(soi_task(coiidx:length(soi_task)-coiidx),'omitnan');
        
        
    end  %end pair loop
    
    % % % % % do sign test window condition % % % % % %
            for j=1:nWinTask
                for p=1:nPairs
                    diff_task_win(j, p)= c_task_win(j,p)-c_rest_tot(p);
                end
            end         
    diffvec=diff_task_win(:);
    [p,~,stats]= signtest(diffvec);
    results(1,win,it)=p; %(p,zval,sign, total);
    results(2,win,it)=stats.zval;
    results(3,win,it)=stats.sign; %number of positive signs
    results(4,win,it)=nWinTask*nPairs; %total number of signs
    
    
end %end window loop


% % % % % do sign test total conditions % % % % % %

%total data, online fashion
diff_task_tot=c_task_tot-c_rest_tot;
diffvec=diff_task_tot;
[p,~,stats]= signtest(diffvec);
results_tot(1,it)=p; %(p,zval,sign, total);
results_tot(2,it)=stats.zval;
results_tot(3,it)=stats.sign; %number of positive signs
results_tot(4,it)=nPairs; %total number of signs


%total data, offline fashion
diff_task_off=c_task_tot_offmat-c_rest_tot_offmat;
diffvec=diff_task_off;
[p,~,stats]= signtest(diffvec);
results_off(1,it)=p; %(p,zval,sign, total);
results_off(2,it)=stats.zval;
results_off(3,it)=stats.sign; %number of positive signs
results_off(4,it)=nPairs; %total number of signs


end %end iterations loop


%% % % % % % make histogram plot % % % % %

ges_results=results;
ges_results(:,nTestWindows+1,:)=results_tot;
ges_results(:,nTestWindows+2,:)=results_off;

ratios=zeros(size(ges_results,2,3));


for i=1:size(ges_results,2)
    for it=1:nIterations
    ratios(i,it)=ges_results(3,i,it)/ges_results(4,i,it);
    end
end
    
mean_ratios=mean(ratios,2);
std_ratios=std(ratios,0,2);
    
for i=1:size(ges_results,2)    
    if i<=nTestWindows
        windows(i)={[num2str(round(windowsizelist(i)/fs)) ' sec']};
    elseif i==nTestWindows+1
        windows(i)={'complete'};
    elseif i==nTestWindows+2
        windows(i)={'offline'};
    end
end
nexttile([2 2]);
X=categorical(windows);
X=reordercats(X,windows); %apparently Matlab is reordering vector when transforming to categorical... need to change back.
bar(X,mean_ratios,'LineWidth', 3); hold on
er=errorbar(X,mean_ratios,std_ratios,std_ratios,'LineWidth', 3);
er.Color=[0 0 0];
er.LineStyle='none';
axes=gca;
axes.FontSize=20;
axes.FontWeight='b';
axes.LineWidth=3;
axes.YTick=[0.25 0.5 0.75 1];
axes.TickLength = [0.05 0.035];
ylim([0 1.1]);
xlabel('window size');

title({'ratio of correctly classified windows'}, 'FontSize', 22, 'FontWeight', 'b');

legend('ratio for 25 sim. pairs', 'std. over 50 iterations','FontSize', 15,'Location','east');



