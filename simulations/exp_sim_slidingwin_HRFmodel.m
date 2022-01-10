%% simulated sliding-window sync of artifical subjects with same model hrf-convolved
%
%KK 13.10.21 adapt for nice figures (for JGM)
%
%KK 31.08.21 script based on HRF_trials_coh.m; add rest condition & sliding 1-step
%window
%
%
% needs SPM hrf function
%


%create signal time course
cond= 'FP' ;
model='peak';
nTrials=30; %number of trials
tRest=60; %60 seconds rest
fs=7.8; %sampling rate 7.8 Hz

switch cond  
    
    case 'FP'
        iti=8; %8sec intertrial interval
        sTask=nTrials*iti*fs; %samples task
        sIti=iti*fs;
        tVec=0:1/fs:iti*nTrials; 
        
    case 'PD/PS'        
        % 60 trials, 60 hrfs (standard SPM hrf), 11sec intertrial interval, fs=7.8
        iti=11;
        sTask=nTrials*iti*fs;
        sIti=iti*fs; %15*7.8=117; mod(0,117)=0, first HR at first time point
        tVec=0:1/fs:iti*nTrials; %for display later
end

tExp=tRest+nTrials*iti;
exp_tVec=0:1/fs:tExp;
sExp=length(exp_tVec);
rest_sampleVec=zeros(tRest*fs,1);
rest_tVec=0:1/fs:tRest;
%sampleVec=zeros(sTask+1,1); %+1 to get matching dims later for display with tVec
modelVec=zeros(sTask+1,1);
HRF=spm_hrf(1/fs);


% %create artifical time course. (feel free to optimize the following code)
% for i=1:sTask
%    if mod(i-1,round(iti))==0 
%        if i-1+length(HRF)<length(sampleVec)
%        dummy=sampleVec(i:i-1+length(HRF));
%        sampleVec(i:i-1+length(HRF))=dummy+HRF;
%        else
%            c=1;
%        for j=i:length(sampleVec)
%            sampleVec(j)=sampleVec(j)+HRF(c);
%            c=c+1;        
%        end
%        end
%    end  
% end
% stvec=(sampleVec-mean(sampleVec))/std(sampleVec);%standardize time course for simplicity

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
        %block model
        blockdur=4; %duration of each block in seconds
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

conv_modelVec=conv(modelVec,HRF);%,'same'); %to get correct convolution, use zero-padded version, later cut off that part
%figure; plot(modelVec); hold on; plot(conv_modelVec); hold off

std_conv_modelVec=(conv_modelVec-mean(conv_modelVec))/std(conv_modelVec); %standardize time course for simplicity

%% create exemplary subjects and get their basic parameters
noiselvl=7;
SNR=20*log(1/noiselvl);

%get params from exemplary subjects

%artifical subject 1
sub1=[randn(size(rest_sampleVec))*noiselvl; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselvl]; %add random gaussian noise
%artifical subject 2
sub2=[randn(size(rest_sampleVec))*noiselvl; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselvl]; %add random gaussian noise


%calculate basic coherence
[Rsq,period,~,coi,~]=wtc([exp_tVec', sub1(1:sExp)],[exp_tVec', sub2(1:sExp)],'mcc',0, 'ArrowSize',0.5, 'ArrowDensity', [30,20]);
for j=1:1:length(coi)
    Rsq(period >= coi(j), j) = NaN; 
end

%get cutoff
poi=[5 11];%[6.6 14.02]; %period of interest in sec
pnoi(1) = find(period > poi(1), 1, 'first');
pnoi(2) = find(period < poi(2), 1, 'last');
%find cutoff (coi) at pnoi
pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
cutoff=find(exp_tVec>pcoi,1,'first');


%% simulate a bunch of additional subjects

nPairs=30;
pairMat=NaN(nPairs,2,length(sub1)); %because of convolution, length(sub1))>sExp; will be cut later

windowsize=30*fs; %30 sec
sWin=sExp-length(rest_sampleVec)-windowsize; %number of data points for sliding window coherence
c_deltaMat=NaN(nPairs,sWin);

pRsq=NaN(nPairs,size(Rsq,1),size(Rsq,2));

for p=1:nPairs
    
    %sub1
    pairMat(p,1,:)=[randn(size(rest_sampleVec))*noiselvl; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselvl];
    %sub2
    pairMat(p,2,:)=[randn(size(rest_sampleVec))*noiselvl; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselvl];
    
    %calc coherence
    [pRsq(p,:,:),period,~,coi,~]=wtc([exp_tVec', squeeze(pairMat(p,1,1:sExp))],[exp_tVec', squeeze(pairMat(p,2,1:sExp))],'mcc',0, 'ArrowSize',0.5, 'ArrowDensity', [30,20]);
    for j=1:1:length(coi)
        pRsq(p,period >= coi(j), j) = NaN; 
    end
    
    %get signal of interest
    soi=squeeze(mean(pRsq(p,pnoi(1):pnoi(2),:),2,'omitnan'));
    
    %calc coherence of rest block
    c_rest=mean(soi(1:length(rest_tVec)),'omitnan');
    
    %calc sliding window coherence
    c_task_win=mean_wtc_sliding_window(soi(tRest*fs+1:end), cutoff, windowsize);
    
    %calc difference task-rest for each window
    c_deltaMat(p,:)=c_task_win-c_rest;

end


c_deltamean=mean(c_deltaMat,1,'omitnan');
c_deltastd=std(c_deltaMat,0,1,'omitnan')/sqrt(size(c_deltaMat,1)); %SEM

Rsqmean=squeeze(mean(pRsq,1,'omitnan'));
Rsqstd=squeeze(std(pRsq,0,1,'omitnan')/sqrt(size(pRsq,1))); %SEM


%% make figures
f=figure;
f.Position(3:4)=[500 800];

tls=tiledlayout(11,2);
tls.TileSpacing='none';
tls.Padding='compact';
title(tls, ['RPS simulation, ' model ' model, sliding window size ' num2str(windowsize/fs) ' sec, SNR=' num2str(SNR) ' dB (=ratio 1/' num2str(noiselvl) ')'], 'FontSize', 16);

%model
ax1=nexttile([3,2]);

total_modelVec=[zeros(length(rest_tVec)-1,1); modelVec(1:end)];
total_conv_modelVec=[zeros(length(rest_tVec)-1,1); std_conv_modelVec(1:length(modelVec))];
ax1.YLim=[-2 1.5];
plot(ax1,exp_tVec, total_modelVec, 'LineWidth',1.5); hold on; plot(ax1,exp_tVec,total_conv_modelVec, 'LineWidth',1.2);hold off;

%plot(ax1,tVec, modelVec); hold on; plot(ax1,tVec,conv_modelVec(1:length(modelVec)));hold off;
%ax1.YLim=[-0.1 1.1];
xlabel('time (sec)');
legend(ax1,{'model','hrf-conv. model'});
ax1.FontSize=13;
title(ax1, 'model');

%artificial subjects
ax2=nexttile([2,2]);
plot(ax2,exp_tVec,sub1(1:length(exp_tVec))); hold on; plot(ax2,exp_tVec,sub2(1:length(exp_tVec))); hold off;
legend(ax2,{'sub1','sub2'});
xlabel('time (sec)');
ax2.FontSize=13;
title(ax2,'two exemplary subjects');

%wtc
ax3=nexttile([3,1]);
wtc([exp_tVec', sub1(1:sExp)],[exp_tVec', sub2(1:sExp)],'mcc',0, 'ArrowSize',0.2, 'ArrowDensity', [30,10])
xlabel('time (sec)');
%wtc(sub1(1:length(exp_tVec)),sub2(1:length(exp_tVec)),'mcc',0, 'ArrowSize',0.2, 'ArrowDensity', [30,10])
%xlabel('time (samples)');
ax3.FontSize=13;
title(ax3, 'exemplary coherence magnitude scalogram');

%frequency/period range
ax5=nexttile([3,1]);
avg=mean(Rsqmean,2,'omitnan');
avgstd=mean(Rsqstd,2,'omitnan'); %for display only, need to find better stat!
periodsec=period;%/fs;
index=find(periodsec<=50,1,'last'); %only interested in periods up to 50 sec (arbitrary)
shadedErrorBar(periodsec(1:index),avg(1:index), avgstd(1:index), 'lineProps',{'-r'});
%plot(ax4,periodsec(1:index),avg(1:index));
%semilogx(periodsec,avg); 
ylabel('mean wtc');
xlabel('period (sec)');
ax5.FontSize=13;
title('average coherence - frequency/period distribution');

%sliding window wtc
ax4=nexttile([3,2]);
c_deltamean_forplot=[zeros(length(rest_sampleVec)+windowsize,1); c_deltamean'];
c_deltastd_forplot=[zeros(length(rest_sampleVec)+windowsize,1); c_deltastd'];
shadedErrorBar(exp_tVec,c_deltamean_forplot, c_deltastd_forplot, 'lineProps',{'-r'});
hold on
hline = refline(0, 0);
hline.Color = 'k';
%slidingwin=[zeros(length(rest_sampleVec)+windowsize,1); c_delta];
%plot(exp_tVec, slidingwin);
ylabel('\Delta wtc');
xlabel('time (sec)');
xline(60);
ax4.YLim=[-0.05 0.15];%0.5];
ax4.FontSize=13;
title(ax4, ['sliding window coherence, window ' num2str(windowsize/fs) ' sec, increment 1 sample, SNR=' num2str(SNR) ' dB (=ratio 1/' num2str(noiselvl) ')']);

linkaxes([ax1,ax2, ax4],'x');
