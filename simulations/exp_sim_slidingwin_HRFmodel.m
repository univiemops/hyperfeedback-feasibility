%% simulated sliding-window sync of artifical subjects with same model hrf-convolved
%
%KK 31.08.21 script based on HRF_trials_coh.m; add rest condition & sliding 1-step
%window, to do: add random additional subjects for group stat
%
%
% needs SPM hrf function
%


%create signal time course
cond= 'PD/PS';
model='block';
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
                    modelVec(i+j)=1;
                end
            end
        end
end

conv_modelVec=conv(modelVec,HRF);%,'same'); %to get correct convolution, use zero-padded version, later cut off that part
%figure; plot(modelVec); hold on; plot(conv_modelVec); hold off

std_conv_modelVec=(conv_modelVec-mean(conv_modelVec))/std(conv_modelVec); %standardize time course for simplicity

%%


noiselevel=2;
%artifical subject 1
sub1=[randn(size(rest_sampleVec))*noiselevel; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselevel]; %add random gaussian noise
%artifical subject 2
sub2=[randn(size(rest_sampleVec))*noiselevel; ...
    std_conv_modelVec+randn(size(std_conv_modelVec))*noiselevel]; %add random gaussian noise






%%

%calculate basic coherence
[Rsq,period,~,coi,~]=wtc([exp_tVec', sub1(1:sExp)],[exp_tVec', sub2(1:sExp)],'mcc',0, 'ArrowSize',0.5, 'ArrowDensity', [30,20]);
for j=1:1:length(coi)
    Rsq(period >= coi(j), j) = NaN; 
end

%get cutoff
poi=[6.6 14.02]; %period of interest in sec
pnoi(1) = find(period > poi(1), 1, 'first');
pnoi(2) = find(period < poi(2), 1, 'last');
%find cutoff (coi) at pnoi
pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
cutoff=find(exp_tVec>pcoi,1,'first');

soi=mean(Rsq(pnoi(1):pnoi(2),:),1,'omitnan'); %signal of interest

c_rest=mean(soi(1:length(rest_tVec)),'omitnan');

windowsize=30*fs; %30 sec
c_task_win=mean_wtc_sliding_window(soi(tRest*fs+1:end), cutoff, windowsize);

c_delta=c_task_win-c_rest;

%% make figures

figure;
f=tiledlayout(5,1);
title(f, ['condition ' cond ', ' model ' model']);

%model
ax1=nexttile;

total_modelVec=[zeros(length(rest_tVec),1); modelVec(2:end)];
total_conv_modelVec=[zeros(length(rest_tVec),1); conv_modelVec(2:length(modelVec))];

plot(ax1,exp_tVec, total_modelVec); hold on; plot(ax1,exp_tVec,total_conv_modelVec);hold off;

%plot(ax1,tVec, modelVec); hold on; plot(ax1,tVec,conv_modelVec(1:length(modelVec)));hold off;
ax1.YLim=[-0.1 1.1];
xlabel('time (sec)');
legend(ax1,{'model','hrf-conv. model'});

%artificial subjects
ax2=nexttile;
plot(ax2,exp_tVec,sub1(1:length(exp_tVec))); hold on; plot(ax2,exp_tVec,sub2(1:length(exp_tVec))); hold off;
legend(ax2,{'sub1','sub2'});
xlabel('time (sec)');


%wtc
ax3=nexttile;
wtc([exp_tVec', sub1(1:sExp)],[exp_tVec', sub2(1:sExp)],'mcc',0, 'ArrowSize',0.2, 'ArrowDensity', [30,10])
xlabel('time (sec)');
%wtc(sub1(1:length(exp_tVec)),sub2(1:length(exp_tVec)),'mcc',0, 'ArrowSize',0.2, 'ArrowDensity', [30,10])
%xlabel('time (samples)');


%sliding window wtc
ax4=nexttile;
slidingwin=[zeros(length(rest_sampleVec)+windowsize,1); c_delta];
plot(exp_tVec, slidingwin);
ylabel('delta wtc');
xlabel('time (sec)');

linkaxes([ax1,ax2,ax3, ax4],'x');


%frequency/period range
ax4=nexttile;
avg=mean(Rsq,2,'omitnan');
periodsec=period/fs;
index=find(periodsec<=50,1,'last'); %only interested in periods up to 50 sec (arbitrary)
plot(ax4,periodsec(1:index),avg(1:index));
%semilogx(periodsec,avg); 
ylabel('mean wtc');
xlabel('period (sec)');