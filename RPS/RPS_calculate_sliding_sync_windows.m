    %%  calculate sync values in different windows
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
    CHOI=15; %channel of interest (if only interested in one channel)


    c=1; %using FP for greatest contrast
    d=2; %number of differences of FP-rest per pair 
    T=30;%trials
    f=7.8125; %Hz, sampling rate
    triallength=floor(8*f); %samples per trial, one trial is 8 sec
    windowlength=triallength*4;%triallength*2; %for now, can change
    
    
  %  Nwd=T/2;%15; %take two trials together 
    


    %coherence
    c_results_rest=zeros(Nch,length(pairs),3)*NaN;
    mc_results_rest=zeros(Nch,length(pairs))*NaN; %for average over all three rest blocks


    c_deltad=zeros(Nch,length(pairs),d,round(T*triallength))*NaN;%for the half trials


    %for coherence - get period of interest poi
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
    
    
    %get passband for filter
    fpass=[1/poi_calc(2) 1/poi_calc(1)];
    
    %find cutoff (coi) at pnoi
    pcoi = max(coi(pnoi(1)),coi(pnoi(2)));
    cutoff=find(t>pcoi,1,'first');
    

    
    
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
            
            
            
            c_results_rest(ch,i,1)=mean(soi(restMat(1,1):restMat(1,2)),'omitnan');
            c_results_rest(ch,i,2)=mean(soi(restMat(2,1):restMat(2,2)),'omitnan');
            c_results_rest(ch,i,3)=mean(soi(restMat(3,1):restMat(3,2)),'omitnan');
            mc_results_rest(ch,i)=mean(c_results_rest(ch,i,:));
 
            

            % and different window lengths, double windows
            xvec=t(subtaskMat(1,1):restMat(2,1))-t(subtaskMat(1,1));
            c_window1=zeros(length(xvec),1)*NaN; %half of mat will actually be empty
            c_window2=c_window1;
            

            


            for v=windowlength:T*triallength%length(xvec)%-windowlength %T*triallength=total number of samples per block
            %first block
            start=subtaskMat(1,1);   
         %   disp(num2str(start));
               assert(start+v<=restMat(2,1), ['sliding window running out of first task block, please check. start+v=' num2str(start+v-1) ', rest=' num2str(restMat(2,1))]);
                    teststat=soi(start+v-windowlength:start+v-1);
            
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end
            c_window1(v)=mean(teststat,'omitnan');  
            
            
            %second block
            start=subtaskMat(1,2);
    %        disp(num2str(start));
               assert(start+v<=restMat(3,1), 'sliding window running out of second task block, please check');
                    teststat=soi(start+v-windowlength:start+v-1);
            for l=1:length(teststat)
                if l<cutoff || l>length(teststat)-cutoff
                    teststat(l)=NaN;
                end
            end
            c_window2(v)=mean(teststat,'omitnan');                




                
                %create delta as measure
                c_deltad(ch,i,1,v)=c_window1(v)-mc_results_rest(ch,i);%c_results_rest(ch,i,1); %first task-first rest
                c_deltad(ch,i,2,v)=c_window2(v)-mc_results_rest(ch,i);%c_results_rest(ch,i,2); %second task-second rest block
                

                end
        end
    end


    
    %% now group stats
    

    
        %c_d=zeros(Nch,length(pairs),Nwd)*NaN;
        %c_d=squeeze(mean(c_deltad,3));

     %for now just plot
     
     %get lines from subtasks
     linevec=squeeze(subtaskMat(:,1))-subtaskMat(1,1);
     
     %get correct x coordinates in sec
     
     
     %xvec=(0:size(c_deltad,4)-1)/fs;
     
     %due to samling slight missmatch, get nearest neighbor in list for
     %plot
     %idx=knnsearch(xvec,linevec');
     
     %%
     
     figure;
     plt=tiledlayout('flow');
     t_window=windowlength/f;
     title(plt,['RPS wtc (task-rest), windowsize ' num2str(t_window) ' sec'], 'FontSize', 16);
     for ch=1:Nch
         c_run1=squeeze(c_deltad(ch,:,1,:));
         c_run1mean=mean(c_run1,1,'omitnan');
         c_run1std=std(c_run1,0,1,'omitnan')/sqrt(size(c_run1,1)); %SEM
         c_run2=squeeze(c_deltad(ch,:,2,:));
         c_run2mean=mean(c_run2,1,'omitnan');
         c_run2std=std(c_run2,0,1,'omitnan')/sqrt(size(c_run2,1));
         
         nexttile
  %          xvec=0:1/f:(length(c_run1mean)-1)/f;%t(subtaskMat(1,1):restMat(2,1))-t(subtaskMat(1,1));%1:length(c_run1mean);

            xvec=1:length(c_run1mean);
         shadedErrorBar(xvec,c_run1mean, c_run1std, 'lineProps',{'-r'}); hold on
         hline = refline(0, 0);
         hline.Color = 'k'; 
         for l=1:length(linevec)
         xline(linevec(l));
         end
         shadedErrorBar(xvec,c_run2mean, c_run2std, 'lineProps',{'-b'}); hold off
         legend('run1', 'run2');
         xlabel('t (sec)')
         ylabel('\Delta wtc')
         title(['channel ' num2str(ch)])
     end
 
   

