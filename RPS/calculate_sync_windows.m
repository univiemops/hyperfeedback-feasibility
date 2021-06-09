    %%  calculate sync values in different windows
    %
    %KK 08.06.21 including pairs with manual markers now
    %
    %KK 19.05.21 debugging
    %
    %KK 04.05.21 implemented correct coi for every step. rest without
    %correct coi, taken from previous analysis.
    %FDR correction also added to double window.
    %basic bandpassfilter for sliding window corr. also added
    %
    %KK 27.04.21 load in wtc for all rest sessions to get as close as
    %possible to paper result. (load in wtc_rest, do not calculate)
    %to do: tidy up code!
    %
    %KK 19.04.21 adding windowed corr, maybe also need to filter?
    %added mean rest
    %to do: figure out good min window length incorp. coi for each window
    %(prob. about 20-30 sec, about 10-15 sec coi miss), code 'dummy coi miss'
    %also at end of window
    %
    %KK 09.04.21 complete code based on own extracted timestamp mat
    %info: use Grinsted's wavelet toolbox
    %take rest as complete block, tasks in steps of subtask
    %compare rest to following task for different steps


    path=['/Users/kathrin/data/RPS/Kathrin_collab/hmrData'];

    pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};
    %pairs_FP_markers=pairs; %having manual markers now
    pairs_FP_markers={'01','02','03','04','05','06','08','09','11','13','14','15','17','18','20','21','22','23','26','28','30','31','32'}; %current workaround, pairs in FP for whom I am actually having markers at the moment. unclear: pair 10, pair 19, pair24, pair 29
    conds={'FP','PS', 'PD','C'};

    Nch=16; %number of channels; channel 1-8 DLPFC, channel 9-16 TPJ
    CHOI=15; %channel of interest (if only interested in one channel)
    Nw=30; %number of windows, for now identical to number of subtasks/trials

    Nwd=15; %take two trials together 

    c=1; %at the moment I am only having markers for FP
    d=2; %number of differences of FP-rest per pair 
    T=30;%triallengths/windowlengths! ramp up
    
    
    %switches
    restload=1; %switch to load in rest data instead of calculating
    filter=0; %switch for filtering for corr
    windows='double';%'single'; %switch for windows
   % artsine='true';
    avgsub=1;
    avgROI=1;

    %correlations
    r_results_rest=zeros(Nch,length(pairs_FP_markers),3)*NaN; %dims: number of channels, number subjects with markers, number of conditions 
    mr_results_rest=zeros(Nch,length(pairs_FP_markers))*NaN; %for average over all three rest blocks
    r_results_task=zeros(Nch,length(pairs_FP_markers),2)*NaN;
    r_results_subtask=zeros(Nch,length(pairs_FP_markers),60)*NaN;
    r_delta=zeros(Nch,length(pairs_FP_markers),d,T)*NaN; %matrix with all delta stats per window length
    r_deltad=zeros(Nch,length(pairs_FP_markers),d,Nwd)*NaN;%for the half trials    
    

    %coherence
    c_results_rest=zeros(Nch,length(pairs_FP_markers),3)*NaN;
    mc_results_rest=zeros(Nch,length(pairs_FP_markers))*NaN; %for average over all three rest blocks
    c_results_task=zeros(Nch,length(pairs_FP_markers),2)*NaN;
    c_results_subtask=zeros(Nch,length(pairs_FP_markers),60)*NaN;
    c_delta=zeros(Nch,length(pairs_FP_markers),d,T)*NaN; %matrix with all delta stats per window length
    c_deltad=zeros(Nch,length(pairs_FP_markers),d,Nwd)*NaN;%for the half trials


    %for coherence - get period of interest poi
    load([path '/Data_' conds{c} '/RPS_01_sub1_' conds{c} '.mat']); %get data of exemplary pair
    hbo1=hbo;
    load([path '/Data_' conds{c} '/RPS_01_sub2_' conds{c} '.mat']); %get data of exemplary pair
    hbo2=hbo;
    poi=[6 20];
    pnoi = zeros(2,1);
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hbo1(:,1)];
    sigPart2 = [t, hbo2(:,1)];
    [~,period,~,coi,~] = wtc(sigPart1, sigPart2, 'mcc', 0);
    pnoi(1) = find(period > poi(1), 1, 'first');
    pnoi(2) = find(period < poi(2), 1, 'last');
    
    poi_calc(1)=period(pnoi(1));
    poi_calc(2)=period(pnoi(2));
    
    
    %get passband for filter
    fpass=[1/poi_calc(2) 1/poi_calc(1)];
    
    %find cutoff (coi) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    
    if restload  %load in rest data for all sessions, don't calc
        wtc_restmat=zeros(length(conds),Nch,length(pairs_FP_markers))*NaN;
%        wtc_restload=zeros(Nch,length(pairs_FP_markers))*NaN;
    end
    

    
    for i=1:length(pairs_FP_markers)
        
        if restload  %load in rest data for all sessions
           for cr=1:length(conds)   
               try
            load(['/Users/kathrin/data/RPS/Kathrin_collab/wtcData/Data_' conds{cr} '/RPS_' pairs_FP_markers{i} '_' conds{cr} '_hbr.mat']);
               catch
               warning(['cant load in /RPS_' pairs_FP_markers{i} '_' conds{cr} '_hbr.mat Leaving NaNs.']);
               continue;
               end
            wtc_restmat(cr,:,i)=coherences(:,2);
           end
           wtc_restload=squeeze(mean(wtc_restmat,1));
        end
        

        load([path '/Data_' conds{c} '/RPS_' pairs_FP_markers{i} '_sub1_' conds{c} '.mat']); %get data (hb and fs)
        hbo1=hbo;
        hbr1=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs_FP_markers{i} '_sub2_' conds{c} '.mat']); %get data (hb and fs)
        hbo2=hbo;
        hbr2=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs_FP_markers{i} '_marker_' conds{c} '.mat']); %get marker matrices

    %add artifical sine wave here
    
        
    %     %correlations
    %     for ch=1:Nch
    %         %for now as proof of principle only corr whole block
    %         %rest 1
    %         r_results_rest(ch,i,1)=corr(hbr1(restMat(1,1):restMat(1,2),ch),hbr2(restMat(1,1):restMat(1,2),ch));
    %         %rest 2
    %         r_results_rest(ch,i,2)=corr(hbr1(restMat(2,1):restMat(2,2),ch),hbr2(restMat(2,1):restMat(2,2),ch));
    %         %rest 3
    %         r_results_rest(ch,i,3)=corr(hbr1(restMat(3,1):restMat(3,2),ch),hbr2(restMat(3,1):restMat(3,2),ch));
    %         %task 1
    %         r_results_task(ch,i,1)=corr(hbr1(taskMat(1,1):taskMat(1,2),ch),hbr2(taskMat(1,1):taskMat(1,2),ch));
    %         %task 2
    %         r_results_task(ch,i,2)=corr(hbr1(taskMat(2,1):taskMat(2,2),ch),hbr2(taskMat(2,1):taskMat(2,2),ch));
    %     end

        %wtc
        %find upper and lower boundaries for period poi
        %get cone of influence coi
        Rsq{Nch} = [];
        Rsq(:) = {NaN(length(period), length(t))};

    for ch=1:Nch %for now do this for every channel separately

            if ~isnan(hbr1(1, ch)) && ~isnan(hbr2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hbr1(:,ch)];
                sigPart2 = [t, hbr2(:,ch)];

                %calc coherence
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence

                %NaN everything outside coi
                for j=1:1:length(coi)
                    Rsq{ch}(period >= coi(j), j) = NaN; %for now use complete calculate, later need coi of substeps
                end
                
                %average over periods to get soi
                soi=mean(Rsq{ch}(pnoi(1):pnoi(2),:),1,'omitnan'); %soi=signal of interest, avaraged over period if interest
%                 
%                 %find cutoff (coi) at pnoi
%                 pcoi = max(coi(pnoi(1)),coi(pnoi(2)));            
                

            end
            
            %now do window calc
            %take rest completely
            
            
            if restload 
                mc_results_rest(ch,i)=wtc_restload(ch,i);
            else
            c_results_rest(ch,i,1)=mean(soi(restMat(1,1):restMat(1,2)),'omitnan');
            c_results_rest(ch,i,2)=mean(soi(restMat(2,1):restMat(2,2)),'omitnan');
            c_results_rest(ch,i,3)=mean(soi(restMat(3,1):restMat(3,2)),'omitnan'); 
            mc_results_rest(ch,i)=mean(c_results_rest(ch,i,:));                
            end
            
            r_results_rest(ch,i,1)=corr(hbr1(restMat(1,1):restMat(1,2),ch),hbr2(restMat(1,1):restMat(1,2),ch), 'rows','complete');
            r_results_rest(ch,i,2)=corr(hbr1(restMat(2,1):restMat(2,2),ch),hbr2(restMat(2,1):restMat(2,2),ch), 'rows','complete');
            r_results_rest(ch,i,3)=corr(hbr1(restMat(3,1):restMat(3,2),ch),hbr2(restMat(3,1):restMat(3,2),ch), 'rows','complete');
            mr_results_rest(ch,i)=mean(r_results_rest(ch,i,:));

            
            % and different window lengths
 
        switch windows
         case 'single'
            c_window1=zeros(Nw,1)*NaN;
            c_window2=c_window1;
            
            r_window1=zeros(Nw,1)*NaN;
            r_window2=r_window1;
            
            for w=1:Nw              
            %do subwindows with correct coi
            teststat=soi(subtaskMat(1,1):subtaskMat(w,1));
%            testhbr1=hbr1(subtaskMat(1,1):subtaskMat(w,1),ch);
%            testhbr2=hbr2(subtaskMat(1,1):subtaskMat(w,1),ch);
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end    
            c_window1(w)=mean(teststat,'omitnan');

            teststat=soi(subtaskMat(1,2):subtaskMat(w,2));
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end             
            c_window2(w)=mean(teststat,'omitnan');

%            c_window1(w)=mean(soi(subtaskMat(1,1):subtaskMat(w,1)),'omitnan');
%            c_window2(w)=mean(soi(subtaskMat(1,2):subtaskMat(w,2)),'omitnan');
            
            r_window1(w)=corr(hbr1(subtaskMat(1,1):subtaskMat(w,1),ch), hbr2(subtaskMat(1,1):subtaskMat(w,1),ch), 'rows','complete');
            r_window2(w)=corr(hbr1(subtaskMat(1,2):subtaskMat(w,2),ch), hbr2(subtaskMat(1,2):subtaskMat(w,2),ch), 'rows','complete');            
    
            %create delta as measure
            c_delta(ch,i,1,w)=c_window1(w)-mc_results_rest(ch,i);%c_results_rest(ch,i,1); %first task-first rest
            c_delta(ch,i,2,w)=c_window2(w)-mc_results_rest(ch,i);%c_results_rest(ch,i,2); %second task-second rest block    

            r_delta(ch,i,1,w)=r_window1(w)-mr_results_rest(ch,i);%r_results_rest(ch,i,1); %first task-first rest %take care, no z-transform at the moment!
            r_delta(ch,i,2,w)=r_window2(w)-mr_results_rest(ch,i);%r_results_rest(ch,i,2); %second task-second rest block                
            end  

        case 'double'
            % and different window lengths, double windows
            c_window1=zeros(Nwd,1)*NaN;
            c_window2=c_window1;

            r_window1=zeros(Nwd,1)*NaN;
            r_window2=r_window1;
            
            for w=1:Nwd
                
            %do subwindows with correct coi
            teststat=soi(subtaskMat(1,1):subtaskMat(2*w,1));
%            testhbr1=hbr1(subtaskMat(1,1):subtaskMat(w,1),ch);
%            testhbr2=hbr2(subtaskMat(1,1):subtaskMat(w,1),ch);
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end    
            c_window1(w)=mean(teststat,'omitnan');

            teststat=soi(subtaskMat(1,2):subtaskMat(2*w,2));
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end             
            c_window2(w)=mean(teststat,'omitnan');                  
                
%             c_window1(w)=mean(soi(subtaskMat(1,1):subtaskMat(2*w,1)),'omitnan');
%             c_window2(w)=mean(soi(subtaskMat(1,2):subtaskMat(2*w,2)),'omitnan');
if filter
    %filter for corr, start with this quick-and-dirty way
    fhbr1 = bandpass(hbr1,fpass,fs); %1 instead of fs
    fhbr2 = bandpass(hbr2,fpass,fs);
    r_window1(w)=corr(fhbr1(subtaskMat(1,1):subtaskMat(2*w,1),ch), fhbr2(subtaskMat(1,1):subtaskMat(2*w,1),ch), 'rows','complete');
    r_window2(w)=corr(fhbr1(subtaskMat(1,2):subtaskMat(2*w,2),ch), fhbr2(subtaskMat(1,2):subtaskMat(2*w,2),ch), 'rows','complete');
else
    r_window1(w)=corr(hbr1(subtaskMat(1,1):subtaskMat(2*w,1),ch), hbr2(subtaskMat(1,1):subtaskMat(2*w,1),ch), 'rows','complete');
    r_window2(w)=corr(hbr1(subtaskMat(1,2):subtaskMat(2*w,2),ch), hbr2(subtaskMat(1,2):subtaskMat(2*w,2),ch), 'rows','complete');
end
            

            %create delta as measure
            c_deltad(ch,i,1,w)=c_window1(w)-mc_results_rest(ch,i);%c_results_rest(ch,i,1); %first task-first rest
            c_deltad(ch,i,2,w)=c_window2(w)-mc_results_rest(ch,i);%c_results_rest(ch,i,2); %second task-second rest block   
            
            r_deltad(ch,i,1,w)=r_window1(w)-mr_results_rest(ch,i);%r_results_rest(ch,i,1); %first task-first rest %take care, no z-transform at the moment!
            r_deltad(ch,i,2,w)=r_window2(w)-mr_results_rest(ch,i);%r_results_rest(ch,i,2); %second task-second rest block
            end   
            
        end

    %         [wcoh,wcs,fcoh]=wcoherence(hbo1(:,ch),hbo2(:,ch),fs); 
    %         c_results_rest(ch,i,1)=mean(wcoh(55,restMat(1,1):restMat(1,2))); %freq 55 as frequency of interest %need to check: cone of influence!
    %         c_results_rest(ch,i,2)=mean(wcoh(55,restMat(2,1):restMat(2,2)));
    %         c_results_rest(ch,i,3)=mean(wcoh(55,restMat(3,1):restMat(3,2)));
    %         c_results_task(ch,i,1)=mean(wcoh(55,taskMat(1,1):taskMat(1,2)));
    %         c_results_task(ch,i,2)=mean(wcoh(55,taskMat(2,1):taskMat(2,2)));
        end


    end

%% now group stats

switch windows
    case 'single'
%- single window
%for now treat both sessions per pair independently (2*N subject pairs)
%->data matrix channels x pairs x windows

c_d=reshape(c_delta,Nch,2*length(pairs_FP_markers),30);
r_d=reshape(r_delta,Nch,2*length(pairs_FP_markers),30);

c_results=zeros(Nch,Nw,2)*NaN;
r_results=c_results;

for i=1:Nch
    for j=1:Nw %wilcoxon signed rank for every window
        if isempty(find(isnan(c_d(i,:,j))==0))  %#ok<*EFIND>
            continue
        end
       [p,h,stats] = signrank(squeeze(c_d(i,:,j)),0,'tail','right');
        c_results(i,j,1)=stats.zval;
        c_results(i,j,2)=p;
        if isempty(find(isnan(r_d(i,:,j))==0)) 
            continue
        end        
       [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
        r_results(i,j,1)=stats.zval;
        r_results(i,j,2)=p;        
    end
end



%create heatmap of p-values
figure;
nexttile;
h=heatmap(squeeze(c_results(:,:,2)));
h.Title='p-values wtc (uncorr., wilcoxon signed rank) for task>rest';
h.XLabel='time windows (integrating), step=8sec';
h.YLabel='channel';
h.Colormap=parula;
h.ColorScaling='log';

nexttile;
h=heatmap(squeeze(r_results(:,:,2)));
h.Title='p-values corr (uncorr., wilcoxon signed rank) for task>rest';
h.XLabel='time windows (integrating), step=8sec';
h.YLabel='channel';
h.Colormap=parula;
h.ColorScaling='log';




    case 'double'
% double length window
%for now treat both sessions per pair independently (2*N subject pairs)
%->data matrix channels x pairs x windows

if ~avgsub

c_d=reshape(c_deltad,Nch,2*length(pairs_FP_markers),Nwd);
r_d=reshape(r_deltad,Nch,2*length(pairs_FP_markers),Nwd);

c_results=zeros(Nch,Nwd,2)*NaN;
r_results=c_results;
for i=1:Nch
    for j=1:Nwd %wilcoxon signed rank for every window
        if isempty(find(isnan(c_d(i,:,j))==0))  %#ok<*EFIND>
            continue
        end
       [p,h,stats] = signrank(squeeze(c_d(i,:,j)),0,'tail','right');
        c_results(i,j,1)=stats.zval;
        c_results(i,j,2)=p;
        if isempty(find(isnan(r_d(i,:,j))==0)) 
            continue
        end        
       [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
        r_results(i,j,1)=stats.zval;
        r_results(i,j,2)=p;        
    end
end

%insert FDR correction here
%A = A(:,~all(isnan(A)));
c_resultsr=squeeze(c_results(:,:,2));
c_resultsr=c_resultsr(:,~all(isnan(c_resultsr)));
c_results_vec=c_resultsr(:);
[n_signif,index_signif]=fdr(c_results_vec,0.05);


%create heatmap of p-values

figure;
nexttile;
h=heatmap(squeeze(c_results(:,:,2)));
h.Title='p-values wtc (uncorr., wilcoxon signed rank) for task>rest';
h.XLabel='time windows (integrating), step=2*8sec';
h.YLabel='channel';
h.Colormap=parula;
h.ColorScaling='log';

%  nexttile;
%  h=heatmap(squeeze(r_results(:,:,2)));
%  h.Title='p-values correlation (uncorr., wilcoxon signed rank) for task>rest';
%  h.XLabel='time windows (integrating), step=2*8sec';
%  h.YLabel='channel';
%  h.Colormap=parula;
%  h.ColorScaling='log';
end






if avgsub
   %c_d=zeros(Nch,length(pairs_FP_markers),Nwd)*NaN;
   c_d=squeeze(mean(c_deltad,3));
for i=1:Nch
    for j=1:Nwd %wilcoxon signed rank for every window
        if isempty(find(isnan(c_d(i,:,j))==0))  %#ok<*EFIND>
            continue
        end
       [p,h,stats] = signrank(squeeze(c_d(i,:,j)),0,'tail','right');
        c_results(i,j,1)=stats.zval;
        c_results(i,j,2)=p;
%         if isempty(find(isnan(r_d(i,:,j))==0)) 
%             continue
%         end        
%        [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
%         r_results(i,j,1)=stats.zval;
%         r_results(i,j,2)=p;        
    end
end

% %insert FDR correction here
% %A = A(:,~all(isnan(A)));
% c_resultsr=squeeze(c_results(:,:,2));
% c_resultsr=c_resultsr(:,~all(isnan(c_resultsr)));
% c_results_vec=c_resultsr(:);
% [n_signif,index_signif]=fdr(c_results_vec,0.05);


%create heatmap of p-values

figure;
nexttile;
h=heatmap(squeeze(c_results(:,:,2)));
h.Title='p-values wtc (uncorr., wilcoxon signed rank) for task>rest';
h.XLabel='time windows (integrating), step=2*8sec';
h.YLabel='channel';
h.Colormap=parula;
h.ColorScaling='log';
   
   
end    

end
