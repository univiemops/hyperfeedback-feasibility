%% time windows with music data
%
%KK 17.08. adapting for sliding wtc windows. throwing out corr for now.
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
%T= 36; %how many windows?  36=about(6600-1800)/fs/(2*phrase) for complete block
%T=32; %31=about(6600-1800-600)/fs/(2*phrase) for block minus learning init
%winsize=round(2*phrase);%roughly 2*phrase=13 sec (multiples of phrase?)
offset=600; 
fs=10; %sampling rate 10 Hz
task_indices=1800+600:6600; %signal of interest for task block
windowsize=300;

defacto_length_results=length(task_indices)-windowsize;
       
    %coherence
    c_results_rest=NaN(Nch,length(pairs));
    c_results_task=NaN(Nch,length(pairs),defacto_length_results);
    c_delta=NaN(Nch,length(pairs),defacto_length_results); %matrix with all delta stats per window length
%    c_deltad=zeros(Nch,length(pairs),d,Nwd)*NaN;%for the half trials




%get exact parameters for period and filter
    poi=[6.6 14.02]; %period of interest

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
                [Rsq{ch}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % Rsq=r square - measure for coherence

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
            
            
            soi_task=soi(task_indices); %signal of interest for task block
            
            
            c_task_win=mean_wtc_sliding_window(soi_task, cutoff, windowsize);
            
            c_results_task(ch,i,:)=c_task_win;

            %create delta as measure
                if ~isfinite(c_task_win)
                    warning(['pair ' num2str(i) 'channel ' num2str(ch) 'needs to be checked.']);
                    c_delta(ch,i,:)=NaN(defacto_length_results,1);
                else
                    c_delta(ch,i,:)=c_results_task(ch,i,:)-c_results_rest(ch,i);
                end
            

            
        end %end of loop over channels
        
        
    end %end of loop over pairs

    
    
%% group stats    

culprits=find(~isfinite(c_delta));
c_delta_noshit=c_delta;
c_delta_noshit(culprits)=NaN;  


c_delta_clean=NaN(size(c_delta_noshit));



c_results_p=NaN(Nch,defacto_length_results);

for c=1:Nch
    for l=1:defacto_length_results %wilcoxon signed rank for every window
        
        vec_without_outliers=rmoutliers(squeeze(c_delta_noshit(c,:,l)),'mean');
        c_delta_clean(c,1:length(vec_without_outliers),l)=vec_without_outliers;
        
        if isempty(find(isnan(c_delta(c,:,l))==0))  %#ok<*EFIND>
            continue
        end
       [p,h,stats] = signrank(squeeze(c_delta(c,:,l)),0,'tail','right');
%        c_results(i,j,1)=stats.zval;
        c_results_p(c,l)=p;
%         if isempty(find(isnan(r_d(i,:,j))==0)) 
%             continue
%         end        
%        [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
%         r_results(i,j,1)=stats.zval;
%         r_results(i,j,2)=p;        
    end
end
    
c_disp=squeeze(mean(c_delta_clean,2));

task_times=(task_indices-(1800+600)+windowsize)/fs;
task_times=task_times(1:end-windowsize);
windowsize_sec=windowsize/fs;

 %%
 figure;
 t=tiledlayout(2,1);
 
 nexttile
 
 clims=[-0.2 0.2];
 Im=imagesc(c_disp,clims);
 ax=gca;
 ax.XTick=1:600:4000;
 strval=mat2str(windowsize/fs:60:4000/fs);
 strval=erase(strval,'[');
 strval=erase(strval,']');
 strval=split(strval,' ');
 ax.XTickLabel=strval;
 ax.FontSize = 15;
 %set(Im, 'XData',task_times);
 c=colorbar;
 c.Label.String = '\Delta wtc';
 xlabel('time (seconds)');
 ylabel('channels');
  title(['wtc(windowed-task)-wtc(total-rest), windowsize=' num2str(windowsize_sec) ' sec']);
 
 nexttile

 clims=[0.001 0.05];
 Im2=imagesc(c_results_p, clims);
  ax=gca;
 ax.XTick=1:600:4000;
 strval=mat2str(windowsize/fs:60:4000/fs);
 strval=erase(strval,'[');
 strval=erase(strval,']');
 strval=split(strval,' ');
 ax.XTickLabel=strval;
 ax.FontSize = 15;
 c=colorbar;
  c.Label.String = 'p (uncorr.)';
  xlabel('time (seconds)');
 ylabel('channels');
  title('p-values (Wilcoxon signed-rank one-sided)');
  
  title(t,'sliding window wtc in incrementing steps of 1 sample', 'FontSize', 18, 'FontWeight', 'bold');
% 
%  h=heatmap(c_disp);
% h.Colormap=parula;
% h.ColorScaling='scaled';   
% h.GridVisible='off';
% h.ColorLimits=[-0.2 0.2];
% 
%  nexttile;
%  h=heatmap(c_results_p);
% h.Colormap=parula;
% h.ColorScaling='scaled';   
% h.GridVisible='off';
% h.ColorLimits=[0.0001 0.05];



%%
% %% old group stats
%  
% c_results=zeros(Nch,T,2)*NaN;
% 
%     
% for i=1:Nch
%     for j=1:T %wilcoxon signed rank for every window
%         if isempty(find(isnan(c_delta(i,:,j))==0))  %#ok<*EFIND>
%             continue
%         end
%        [p,h,stats] = signrank(squeeze(c_delta(i,:,j)),0,'tail','right');
% %        c_results(i,j,1)=stats.zval;
%         c_results(i,j,2)=p;
% %         if isempty(find(isnan(r_d(i,:,j))==0)) 
% %             continue
% %         end        
% %        [p,h,stats] = signrank(squeeze(r_d(i,:,j)),0,'tail','right');
% %         r_results(i,j,1)=stats.zval;
% %         r_results(i,j,2)=p;        
%     end
% end
% 
% % %insert FDR correction here
% %A = A(:,~all(isnan(A)));
% c_resultsr=squeeze(c_results(:,:,2));
% c_resultsr=c_resultsr(:,~all(isnan(c_resultsr)));
% c_results_vec=c_resultsr(:);
% [n_signif,index_signif]=fdr(c_results_vec,0.05);
% c_fdrvec=zeros(size(c_results_vec))*NaN;
% c_fdrvec(index_signif)=c_results_vec(index_signif);
% c_fdr=reshape(c_fdrvec, size(c_resultsr));
% 
% %c_fdr_dirtydisplay=[NaN(16,1) c_fdr];
% 
% 
% %create heatmap of p-values
% 
% figure;
% nexttile;
% h=heatmap(squeeze(c_results(:,:,2)));
% h.Title='p-values wtc (uncorr., wilcoxon signed rank) for task>rest';
% h.XLabel='time windows (integrating), step=2*6.5sec';
% h.YLabel='channel';
% h.Colormap=parula;
% h.ColorScaling='scaled';%'log';
% h.ColorLimits=[0.0001 0.05];
% h.FontSize=16;
%              
%    %heatmap only for FDR corr.
% figure;
% nexttile;
% h=heatmap(c_fdr);
% h.Title='p-values wtc (uncorr., wilcoxon signed rank) for task>rest';
% h.XLabel='time windows (integrating), step=2*6.5sec';
% h.YLabel='channel';
% h.Colormap=parula;
% h.ColorScaling='scaled';%'log';
% h.ColorLimits=[0.0001 0.05];
% h.FontSize=16;         
% % 
% %  