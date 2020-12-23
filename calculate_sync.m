%calculate sync values in different blocks

%get hmr values plus markers
%calculate sync from start marker to end marker
%get values per subject
%group stats?

path=['/Users/kathrin/data/RPS/Kathrin_collab/hmrData'];

pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};
markers_FP={'01','02','03','04','05','06','08','09','11','13','14','15','17','18','19','20','21','22','23','26','28','30','31','32'}; %current workaround, pairs in FP for whom I am actually having markers at the moment.
conds={'FP','PS', 'PD','C'};

Nch=16; %number of channels; channel 1-8 DLPFC, channel 9-16 TPJ
CHOI=15; %channel of interest (if only interested in one channel)

c=1; %at the moment I am only having markers for FP

r_results_rest=zeros(Nch,length(markers_FP),3)*NaN;
r_results_task=zeros(Nch,length(markers_FP),2)*NaN;

c_results_rest=zeros(Nch,length(markers_FP),3)*NaN;
c_results_task=zeros(Nch,length(markers_FP),2)*NaN;

%for dummy stats
r_rest=zeros(length(markers_FP),1)*NaN;
r_task=zeros(length(markers_FP),1)*NaN;

c_rest=zeros(length(markers_FP),1)*NaN;
c_task=zeros(length(markers_FP),1)*NaN;

for i=1:length(markers_FP)
    
    load([path '/Data_' conds{c} '/RPS_' markers_FP{i} '_sub1_' conds{c} '.mat']); %get data (hb and fs)
    hbo1=hbo;
    hbr1=hbr;
    load([path '/Data_' conds{c} '/RPS_' markers_FP{i} '_sub2_' conds{c} '.mat']); %get data (hb and fs)
    hbo2=hbo;
    hbr2=hbr;
    load([path '/Data_' conds{c} '/RPS_' markers_FP{i} '_marker_' conds{c} '.mat']); %get marker matrices
     
    for ch=1:Nch
        %for now as proof of principle only corr whole block
        %rest 1
        r_results_rest(ch,i,1)=corr(hbr1(restMat(1,1):restMat(1,2),ch),hbr2(restMat(1,1):restMat(1,2),ch));
        %rest 2
        r_results_rest(ch,i,2)=corr(hbr1(restMat(2,1):restMat(2,2),ch),hbr2(restMat(2,1):restMat(2,2),ch));
        %rest 3
        r_results_rest(ch,i,3)=corr(hbr1(restMat(3,1):restMat(3,2),ch),hbr2(restMat(3,1):restMat(3,2),ch));
        %task 1
        r_results_task(ch,i,1)=corr(hbr1(taskMat(1,1):taskMat(1,2),ch),hbr2(taskMat(1,1):taskMat(1,2),ch));
        %task 2
        r_results_task(ch,i,2)=corr(hbr1(taskMat(2,1):taskMat(2,2),ch),hbr2(taskMat(2,1):taskMat(2,2),ch));
        
        %wtc
        [wcoh,wcs,fcoh]=wcoherence(hbo1(:,ch),hbo2(:,ch),fs); 
        c_results_rest(ch,i,1)=mean(wcoh(55,restMat(1,1):restMat(1,2))); %freq 55 as frequency of interest %need to check: cone of influence!
        c_results_rest(ch,i,2)=mean(wcoh(55,restMat(2,1):restMat(2,2)));
        c_results_rest(ch,i,3)=mean(wcoh(55,restMat(3,1):restMat(3,2)));
        c_results_task(ch,i,1)=mean(wcoh(55,taskMat(1,1):taskMat(1,2)));
        c_results_task(ch,i,2)=mean(wcoh(55,taskMat(2,1):taskMat(2,2)));
    end
    
    %dummy stats just as proof of principle (use only TPJ channels)
    r_rest(i)=mean(mean(r_results_rest(9:16,i,:)));
    r_task(i)=mean(mean(r_results_task(9:16,i,:)));
    
    c_rest(i)=mean(mean(c_results_rest(9:16,i,:)));
    c_task(i)=mean(mean(c_results_task(9:16,i,:)));    
end