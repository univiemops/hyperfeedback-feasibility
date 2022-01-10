    %%  calculate sync values in different windows
    %
    %KK 08.01. calc CNR instead of delta wtc
    % note that there is a maximum resolution for calculating CNR
    %since a rest block is only ~65 sec and variability of the window calc
    %in the rest block has to be evaluated. Taking max. ~55 sec (7 trials).
    %taking mean of variability in all three rest blocks for now
    %to do: final figure, group stuff
    %
    %
    %KK 15.12. loop over windowsizes
    %
    %KK 13.10. make nice figures for JGM
    %
    %KK 27.07. improve plot. Need to do: upload debugged get-markers scrpt
    %
    %
    %KK 26.07. create a real sliding window increasing in steps of one
    %sample, debug
    %
    %
    %KK 21.07. make new slimmer script just for wtc windows of double size.
    %removed: restload, windowswitch, assumed avgsub
    %take different starting points
    %
    %KK 30.06. to do: z-score corr already during averaging of rest! check filter for
    %rest!!!!!
    %
    %KK 24.06.21 alternative frequency range: [5 11] - before [6 20]. BUT: would need to
    %recalculate wtc for all rest conditions!
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

    conds={'FP','PS', 'PD','C'};

    Nch=16; %number of channels; channel 1-8 DLPFC, channel 9-16 TPJ


    c=1; %using FP for greatest contrast
    d=2; %number of differences of FP-rest per pair 
    T=30;%trials
    f=7.8125; %Hz, sampling rate
    triallength=floor(8*f); %samples per trial, one trial is 8 sec
    RunSamp=8*30*f; %exact number
    RestSamp=509; %65.152 sec rest duration ideally
    %windowlength=triallength*8;%triallength*2; %for now, can change
    
    windowsizelist=[triallength*2 triallength*4 triallength*6];
    channellist=[13,14,15,16];
    
    %% loop over windows
    
    for win=1:length(windowsizelist)
    
    windowsize=windowsizelist(win);
    
    defacto_length_results=RunSamp-windowsize;
    defacto_length_rest=RestSamp-windowsize;
    
    %coherence
 %   c_results_rest=NaN(Nch,length(pairs));
    c_results_rest_win1=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_rest_win2=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_rest_win3=NaN(Nch,length(pairs), defacto_length_rest);
    c_results_task1=NaN(Nch,length(pairs),defacto_length_results);
    c_results_task2=NaN(Nch,length(pairs),defacto_length_results);
    c_delta=NaN(Nch,length(pairs),defacto_length_results*2); %matrix with all delta stats per window length
    

    %for average values
    c_results_rest=zeros(Nch,length(pairs),3)*NaN; 
    mc_results_rest=zeros(Nch,length(pairs))*NaN; %for average over all three rest blocks
    std_avg_rest=NaN(Nch,length(pairs));

    
    c_deltad=zeros(Nch,length(pairs),d,defacto_length_results)*NaN;

    
    %% for coherence - get period of interest poi
    load([path '/Data_' conds{c} '/RPS_01_sub1_' conds{c} '.mat']); %get data of exemplary pair
    hbo1=hbo;
    load([path '/Data_' conds{c} '/RPS_01_sub2_' conds{c} '.mat']); %get data of exemplary pair: hbo, hbr, t (time stamps of measurement)
    hbo2=hbo;
    poi=[5 11];%[6 10];%[5 11]; %[6 20];
    pnoi = zeros(2,1);
    poi_calc=zeros(2,1);
    
    sigPart1 = [t, hbo1(:,1)];
    sigPart2 = [t, hbo2(:,1)];
    [~,period,~,coi,~] = wtc(sigPart1, sigPart2, 'mcc', 0);
    pnoi(1) = find(period > poi(1), 1, 'first');
    pnoi(2) = find(period < poi(2), 1, 'last');
    
    poi_calc(1)=period(pnoi(1));
    poi_calc(2)=period(pnoi(2));
    
    
    
    %find cutoff (coi) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    

    %% loop over dyads
    
    for i=1:length(pairs)
        
        
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub1_' conds{c} '.mat']); %get data (hb,t)
        hbo1=hbo;
        hbr1=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub2_' conds{c} '.mat']); %get data (hb,t)
        hbo2=hbo;
        hbr2=hbr;
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_marker_' conds{c} '.mat']); %get marker matrices
        
        if ~exist('t') || ~exist('fs') ||  ~exist('restMat') ||  ~exist('taskMat') || ~exist('subtaskMat')
            error(['failed to load in markers and/or data for pair ' num2str(pairs{i})]);
        end
        
        %wtc
        %find upper and lower boundaries for period poi
        %get cone of influence coi
        Rsq{Nch} = [];
        Rsq(:) = {NaN(length(period), length(t))};
        
        %% loop over channels of interest
        for chcount=1:length(channellist) %Nch %for now do this for every channel separately
            
            ch=channellist(chcount);
            
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
                
            end
            
            
            %now do window calc
            
            %take rest completely for average value
            c_results_rest(ch,i,1)=mean(soi(restMat(1,1):restMat(1,2)),'omitnan');
            c_results_rest(ch,i,2)=mean(soi(restMat(2,1):restMat(2,2)),'omitnan');
            c_results_rest(ch,i,3)=mean(soi(restMat(3,1):restMat(3,2)),'omitnan');
            mc_results_rest(ch,i)=mean(c_results_rest(ch,i,:));
            
            
            soi_task1=soi(subtaskMat(1,1):restMat(2,1)-1);
            soi_task2=soi(subtaskMat(1,2):restMat(3,1)-1);
            soi_rest1=soi(restMat(1,1):restMat(1,2)-1);
            soi_rest2=soi(restMat(2,1):restMat(2,2)-1);
            soi_rest3=soi(restMat(3,1):restMat(3,2)-1);
            
            
         %   % %
%                    mean_wtc_window_vec=NaN(length(soi_task1)-windowsize,1);
% 
%                     for w=1:length(soi_task1)-windowsize
% 
%                         teststat1=soi_task1(w:w+windowsize-1);
% 
%                         for l=1:length(teststat1)
%                             if l<cutoff || l>length(teststat1)-cutoff
%                                 teststat1(l)=NaN;
%                             end
%                         end
%                         mean_wtc_window_vec(w)=mean(teststat1,'omitnan'); 
%                     end
%                     
%                     
            % %
            
            c_task_win1=mean_wtc_sliding_window(soi_task1, cutoff, windowsize);
            c_task_win2=mean_wtc_sliding_window(soi_task2, cutoff, windowsize); 
            
            c_rest_win1=mean_wtc_sliding_window(soi_rest1, cutoff, windowsize); 
            c_rest_win2=mean_wtc_sliding_window(soi_rest2, cutoff, windowsize);   
            c_rest_win3=mean_wtc_sliding_window(soi_rest3, cutoff, windowsize);
            
            std_avg_rest(ch,i)=(std(c_rest_min1)+std(c_rest_min2)+std(c_rest_min3))/3; %since all rest blocks are of equal length, this is the same as mean adjusting each rest block's windows and then calculating total std
                      
            
            c_deltad(ch,i,1,:)=c_task_win1(1:defacto_length_results)-mc_results_rest(ch,i);  %here are sometimes 2-3 samples missmatch between trigger logging and designed length. Taking designed length to be identical for all sub. For now shouldn't change results but check later!
            c_deltad(ch,i,2,:)=c_task_win2(1:defacto_length_results)-mc_results_rest(ch,i);
            

        end %end channel loop
    end %end dyad loop
    
%% calc CNR

    Cohens_d_dyad=NaN(Nch, length(pairs),2, defacto_length_results);
    
    for c=1:Nch
        for d=1:length(pairs)
            
%            std_rest_avg=(std(c_rest_min1)+std(c_rest_min2)+std(c_rest_min3))/3; %since all rest blocks are of equal length, this is the same as mean adjusting each rest block's windows and then calculating total std
            
            for l=1:defacto_length_results
                for run=1:2
                    
                    %corr for every window for effect size
                    %     c_corr(c,l)=corr(squeeze(c_task_noshit(c,:,l))', squeeze(c_results_rest(c,:))');
                    
                    Cohens_d_dyad(c,d,run,l)= c_delta(c,d,run,l)/std_avg_rest(c,d);
                    
                    %    Cohens_d_std(c,l)=std(c_deltaMat,0,1,'omitnan')/sqrt(size(c_deltaMat,1)); %SEM
                    %     Cohens_d_corrected(c,l)=Cohens_d_simple(c,l)*1/sqrt(2*(1-c_corr(c,l)));
                    
                    %         %wilcoxon signed rank for every window
                    %
                    %         if isempty(find(isnan(c_delta(c,:,l))==0))  %#ok<*EFIND>
                    %             continue
                    %         end
                    %        [p,h,stats] = signrank(squeeze(c_delta(c,:,l)),0,'tail','right');
                    %
                    %         c_results_p(c,l)=p;
                    
                    
                    %clean data for better display
                    %      c_corr_vec_without_outliers=rmoutliers(squeeze(c_corr(c,:,l)),'mean');
                    c_delta_vec_without_outliers=rmoutliers(squeeze(c_delta_noshit(c,:,l)),'mean');
                    c_delta_clean(c,1:length(c_delta_vec_without_outliers),l)=c_delta_vec_without_outliers;
                    
                end
            end
            
        end
    end
    
    
    
    
%end window loop
    
     %% plot stuff
     
     linevec=squeeze(subtaskMat(:,1))-subtaskMat(1,1);
     

     

     fig=figure;
fig.Position(3:4)=[1200 900];
tls=tiledlayout(4,4);
tls.TileSpacing='none';
tls.Padding='compact';

     t_window=windowsize/f;
     title(tls,['RPS wtc (task-rest), windowsize ' num2str(t_window) ' sec'], 'FontSize', 16);
     for ch=1:Nch
         c_run1=squeeze(c_deltad(ch,:,1,:));
         c_run1mean=mean(c_run1,1,'omitnan');
         c_run1std=std(c_run1,0,1,'omitnan')/sqrt(size(c_run1,1)); %SEM
         c_run2=squeeze(c_deltad(ch,:,2,:));
         c_run2mean=mean(c_run2,1,'omitnan');
         c_run2std=std(c_run2,0,1,'omitnan')/sqrt(size(c_run2,1));
         
         ax(ch)=nexttile;
  %          xvec=0:1/f:(length(c_run1mean)-1)/f;%t(subtaskMat(1,1):restMat(2,1))-t(subtaskMat(1,1));%1:length(c_run1mean);

            xvec=1:length(c_run1mean);
            tvectmp=xvec/f;
            tvec=zeros(1,length(xvec));
            tvec(2:end)=tvectmp(1:end-1);
         shadedErrorBar(tvec,c_run1mean, c_run1std, 'lineProps',{'-r'}); hold on
         hline = refline(0, 0);
         hline.Color = 'k'; 
         for l=1:length(linevec)
         xline(linevec(l)/f);
         end
         shadedErrorBar(tvec,c_run2mean, c_run2std, 'lineProps',{'-b'}); hold off
         legend('run1', 'run2');
         xlabel('t (sec)')
         ylabel('\Delta wtc')
         title(['channel ' num2str(ch)], 'FontSize', 15)
         ax(ch).FontSize=13;
         ax(ch).XLim=[0 250];
         
       
     end
     
     linkaxes([ax(1),ax(2),ax(3),ax(4),ax(5),ax(6),ax(7),ax(8),ax(9),ax(10),ax(11),ax(12),ax(13),ax(14),ax(15),ax(16)],'xy');
 
    end
