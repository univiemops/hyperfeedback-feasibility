%% time windows with music data
%
%
%KK 29.06., 30.06. adapting from calculate_sync_windows for RPS
%beware: check faulty channels from preprocessing!
% to do: filter for corr
%
%
%FF is the prefix for the “part learning” condition
%FF2 is the prefix for the “whole learning” condition
%Probe1 = right hemisphere sub 1
%Probe2 = left hem sub 1
%Probe3 = right hem sub 2
%Probe4 = left hem sub2
%time points 1:1800 for resting state and 1801:6600 for task phase
%from paper: 
%FOI: 0.07-0.15Hz  (period 6.61-14.01s)
%musical phrase 6.64\pm 1.56s
% compute complete rest block, take blocks of task phase
% taking hbo, as in paper


path=['/Users/kathrin/data/Music_fNIRS/data'];
pairs={'01','02','03','04','05','06', '07','08','09','10','11','12'};
Nch=44; %bad channels will have zeros
phrase=6.61;%sec
T= 36; %how many windows?  36=about(6600-1800)/fs/(2*phrase)
winsize=round(2*phrase);%roughly 2*phrase=13 sec (multiples of phrase?)
filter=0; %filter for corr

%define vars
    %correlations
    r_results_rest=zeros(Nch,length(pairs))*NaN; %dims: number of channels, number subjects with markers, number of conditions 
    r_delta=zeros(Nch,length(pairs),T)*NaN; %matrix with all delta stats per window length
%    r_deltad=zeros(Nch,length(pairs),Nwd)*NaN;%for the half trials    
    

    %coherence
    c_results_rest=zeros(Nch,length(pairs))*NaN;
    c_delta=zeros(Nch,length(pairs),T)*NaN; %matrix with all delta stats per window length
%    c_deltad=zeros(Nch,length(pairs),d,Nwd)*NaN;%for the half trials




%get exact parameters for period and filter
    poi=[6.6 14.02];

    load([path '/FF_3_MES_Probe1.mat']); %get data of exemplary pair
    hbo1test=hbo;
    load([path '/FF_3_MES_Probe3.mat']); %get data of exemplary pair
    hbo2test=hbo;

    pnoi = zeros(2,1);
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hbo1test(:,14)];
    sigPart2 = [t, hbo2test(:,14)];
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


%    wtc(sigPart1, sigPart2, 'mcc', 0); %have a look at the data

%now real stuff
    for i=1:length(pairs)
        
     %load in participants
        hbo1=zeros(size(hbo1test,1),Nch);
        hbr1=hbo1;
        hbo2=zeros(size(hbo2test,1),Nch);
        hbr2=hbo2;
        %for a start only part learning
        %right hem both participants
        load([path '/FF_' num2str(i) '_MES_Probe1.mat']); %get data (hb and fs)
        hbo1(:,1:Nch/2)=hbo; %first 22 channels right hem
        hbr1(:,1:Nch/2)=hbr;
        load([path '/FF_' num2str(i) '_MES_Probe3.mat']); %get data (hb and fs)
        hhbo2(:,1:Nch/2)=hbo;
        hbr2(:,1:Nch/2)=hbr;        

        %left hem both participants
        load([path '/FF_' num2str(i) '_MES_Probe2.mat']); %get data (hb and fs)
        hbo1(:,Nch/2+1:Nch)=hbo; %next 22 channels left hem
        hbr1(:,Nch/2+1:Nch)=hbr;
        load([path '/FF_' num2str(i) '_MES_Probe4.mat']); %get data (hb and fs)
        hhbo2(:,Nch/2+1:Nch)=hbo;
        hbr2(:,Nch/2+1:Nch)=hbr;     
        
        
        parfor ch=1:Nch %for now do this for every channel separately

            %get signal of interest
            if ~isnan(hbo1(1, ch)) && ~isnan(hbo2(1, ch))   % check if this channel was not rejected in both subjects during preprocessing
                sigPart1 = [t, hbr1(:,ch)];
                sigPart2 = [t, hbr2(:,ch)];

                %calc coherence
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence

                %NaN everything outside coi
                for j=1:1:length(coi)
                    Rsq{ch}(period >= coi(j), j) = NaN; %outer bracket, but for real time need to NaN every window at borders
                end

                %average over periods to get soi
                soi=mean(Rsq{ch}(pnoi(1):pnoi(2),:),1,'omitnan'); %soi=signal of interest, avaraged over period if interest                        
            end
            
            %now do window calc
            %calculate rest block completely
            c_results_rest(ch,i)=mean(soi(1:1800),'omitnan');
            %and different window lengths separately
            c_window=zeros(T,1)*NaN;
            
            
            if filter
            fhbo1 = bandpass(hbo1,fpass,fs); %1 instead of fs
            fhbo2 = bandpass(hbo2,fpass,fs);     
            r_results_rest(ch,i)=corr(fhbo1(1:1800,ch),fhbo2(1:1800,ch),'rows','complete');
            else    
            r_results_rest(ch,i)=corr(hbo1(1:1800,ch),hbo2(1:1800,ch),'rows','complete');
            end
            r_window=zeros(T,1)*NaN;
            
 

            
            for w=1:T
                
            %do subwindows with correct coi
            teststat=soi(1800:1800+w*winsize*fs);
%            testhbr1=hbr1(subtaskMat(1,1):subtaskMat(w,1),ch);
%            testhbr2=hbr2(subtaskMat(1,1):subtaskMat(w,1),ch);
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end    
            c_window(w)=mean(teststat,'omitnan');
           

            if filter
            %filter for corr, start with this quick-and-dirty way
            r_window(w)=corr(fhbo1(1800:1800+w*winsize*fs,ch), fhbo2(1800:1800+w*winsize*fs,ch), 'rows','complete');
            else
            r_window(w)=corr(hbo1(1800:1800+w*winsize*fs,ch), hbr2(1800:1800+w*winsize*fs,ch), 'rows','complete');
            end
            

            %create delta as measure
            if ~isfinite(c_window(w))
                warning(['pair ' num2str(i) 'channel ' num2str(ch) 'needs to be checked for window ' num2str(w)]);
                c_delta(ch,i,w)=NaN;
            else
                c_delta(ch,i,w)=c_window(w)-c_results_rest(ch,i);
            end
            r_delta(ch,i,w)=r_window(w)-r_results_rest(ch,i); %potentially need to z-score but for following wilcoxon no need
            end   
            
        end
        
        
    end     
        
    
%% now group stats
 
c_results=zeros(Nch,T,2)*NaN;
r_results=c_results;
    
for i=1:Nch
    for j=1:T %wilcoxon signed rank for every window
        if isempty(find(isnan(c_delta(i,:,j))==0))  %#ok<*EFIND>
            continue
        end
       [p,h,stats] = signrank(squeeze(c_delta(i,:,j)),0,'tail','right');
%        c_results(i,j,1)=stats.zval;
        c_results(i,j,2)=p;
%         if isempty(find(isnan(r_d(i,:,j))==0)) 
%             continue
%         end        
%        [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
%         r_results(i,j,1)=stats.zval;
%         r_results(i,j,2)=p;        
    end
end

% % %insert FDR correction here
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
h.XLabel='time windows (integrating), step=2*6.5sec';
h.YLabel='channel';
h.Colormap=parula;
h.ColorScaling='log';
             
            

 