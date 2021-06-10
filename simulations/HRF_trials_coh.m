%% try out sync due to similar hrf shape
%
%KK 10.6.21 created artifical timecourse of hrf (assuming delta-shaped
%activation) plus noise. To do: frequency plot not correct yet.
%
%
% needs SPM hrf function
% 60 trials, 60 hrfs (standard SPM hrf), 15sec intertrial interval, fs=7.8
%
fs=7.8;
%t=zeros(4*60*fs,1);
T=60*15*7.8;
tvec=zeros(T,1);
HRF=spm_hrf(1/fs);

%create artifical time course. (feel free to optimize the following code)
for i=1:T
   if mod(i-1,117)==0 %15*7.8=117; mod(0,117)=0, first HR at first time point
       if i-1+length(HRF)<length(tvec)
       dummy=tvec(i:i-1+length(HRF));
       tvec(i:i-1+length(HRF))=dummy+HRF;
       else
           c=1;
       for j=i:length(tvec)
           tvec(j)=tvec(j)+HRF(c);
           c=c+1;        
       end
       end
   end
    
end

stvec=(tvec-mean(tvec))/std(tvec);%standardize time course for simplicity

%artifical subject 1
sub1=stvec+randn(size(stvec)); %add random gaussian noise
%artifical subject 2
sub2=stvec+randn(size(stvec)); %same as sub1 add random gaussian noise

%calculate coherence

[Rsq,period,~,coi,~]=wtc(sub1,sub2,'mcc',0, 'ArrowSize',0.5, 'ArrowDensity', [30,20]);
for j=1:1:length(coi)
    Rsq(period >= coi(j), j) = NaN; 
end

%dimensions are not correct yet.
avg=mean(Rsq,2,'omitnan');
figure; plot(avg);