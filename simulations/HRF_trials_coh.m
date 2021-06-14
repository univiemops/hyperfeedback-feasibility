%% try out sync of two artifical subjects with same model hrf-convolved
%
%KK 14.06.21 generalize to different experimental conditions; can choose
%between delta peaks and blocks
%
%KK 10.06.21 created artifical timecourse of hrf (assuming delta-shaped
%activation) plus gaussian noise. 
%
%
% needs SPM hrf function
%

cond= 'PD/PS';
model='block';

switch cond  
    
    case 'FP'
        fs=7.8;
        iti=8; %8sec
        T=60*iti*fs;
        Iti=iti*fs;
        t=0:1/fs:iti*60; 
        
    case 'PD/PS'        
        % 60 trials, 60 hrfs (standard SPM hrf), 11sec intertrial interval, fs=7.8
        fs=7.8;
        iti=11;
        T=60*iti*fs;
        Iti=iti*fs; %15*7.8=117; mod(0,117)=0, first HR at first time point
        t=0:1/fs:iti*60; %for display later
end

tvec=zeros(T+1,1); %+1 to get matching dims later for display with t
mvec=zeros(T+1,1);
HRF=spm_hrf(1/fs);


% %create artifical time course. (feel free to optimize the following code)
% for i=1:T
%    if mod(i-1,round(iti))==0 
%        if i-1+length(HRF)<length(tvec)
%        dummy=tvec(i:i-1+length(HRF));
%        tvec(i:i-1+length(HRF))=dummy+HRF;
%        else
%            c=1;
%        for j=i:length(tvec)
%            tvec(j)=tvec(j)+HRF(c);
%            c=c+1;        
%        end
%        end
%    end  
% end
% stvec=(tvec-mean(tvec))/std(tvec);%standardize time course for simplicity

%create artifical time course via convolution
switch model
    case'peak'
        %delta peak
        for i=1:T
            if mod(i-1,round(Iti))==0 %15*7.8=117; mod(0,117)=0, first HR at first time point
                mvec(i)=1;
            end
        end
    case 'block'
        %block model
        blockdur=4; %duration of each block in seconds
        boxlength=blockdur*fs;
        for i=1:T
            if mod(i-1,round(Iti))==0 
                for j=1:boxlength
                    mvec(i+j)=1;
                end
            end
        end
end

mhvec=conv(mvec,HRF);%,'same'); %to get correct convolution, use zero-padded version, later cut off that part
%figure; plot(mvec); hold on; plot(mhvec); hold off

stmhvec=(mhvec-mean(mhvec))/std(mhvec); %standardize time course for simplicity

%%
noiselevel=8;
%artifical subject 1
sub1=stmhvec+randn(size(stmhvec))*noiselevel; %add random gaussian noise
%artifical subject 2
sub2=stmhvec+randn(size(stmhvec))*noiselevel; %add random gaussian noise


%calculate coherence
[Rsq,period,~,coi,~]=wtc(sub1,sub2,'mcc',0, 'ArrowSize',0.5, 'ArrowDensity', [30,20]);
for j=1:1:length(coi)
    Rsq(period >= coi(j), j) = NaN; 
end
%% make figures

figure;
f=tiledlayout(4,1);
title(f, ['condition ' cond ', ' model ' model']);

%model
ax1=nexttile;
plot(ax1,t, mvec); hold on; plot(ax1,t,mhvec(1:length(mvec)));hold off;
ax1.YLim=[-0.1 1.1];
xlabel('time (sec)');
legend(ax1,{'model','hrf-conv. model'});

%artificial subjects
ax2=nexttile;
plot(ax2,t,sub1(1:length(tvec))); hold on; plot(ax2,t,sub2(1:length(tvec))); hold off;
legend(ax2,{'sub1','sub2'});
xlabel('time (sec)');
linkaxes([ax1,ax2],'x');

%wtc
ax3=nexttile;
wtc(sub1(1:length(tvec)),sub2(1:length(tvec)),'mcc',0, 'ArrowSize',0.2, 'ArrowDensity', [30,10])
xlabel('time (samples)');

%frequency/period range
ax4=nexttile;
avg=mean(Rsq,2,'omitnan');
periodsec=period/fs;
index=find(periodsec<=50,1,'last'); %only interested in periods up to 50 sec (arbitrary)
plot(ax4,periodsec(1:index),avg(1:index));
%semilogx(periodsec,avg); 
ylabel('mean wtc');
xlabel('period (sec)');
