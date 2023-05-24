%% get markers from files 's' for 'rest' and 'task' (overall task) to mark
%exact starting/ endpoints

%only for condition FP valid

%assume three rest blocks and two task blocks per file

%assume furthermore: (valid at least for FP)
%s(:,1) is rest
%s(:,2) is task
%s(:,3) is also task, more specific task. Called 'subtask' furtheronwards.
%s always demarks starting point of incidents.



%for rest take beginning of rest to task marker
%for task take first and last subtask marker


%loop over all subjects and conditions
%load in 's'
%calculate diffs
%write into matrices to save later on





path=['yourPathToData/hmrData']; %Unix convention used, adapt if on Windows system

pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};

conds={'FP','PS', 'PD','C'};
flag_nRestgrtrnTask=false;

nomatch=zeros(length(pairs),length(conds))*NaN;



    c=1; %use only Free Play 'FP'
    
    errorcount=0;
    errorcount_correctRandT=0;
    errorcount_correctRnotT=0;
    N_done=0;
    
    for i=1:length(pairs) 
        flag_nRestgrtrnTask=false;
        flag_subj2needed=false;
        flag_hopelesscase=false;
        
        restMat=zeros(3,2)*NaN;
        taskMat=zeros(2,2)*NaN;
        subtaskMat=zeros(30,2)*NaN;

   %load sub 1     
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub1_' conds{c} '.mat']);
        %vecs to get all the markers for each condition
        s1_rest=find(s(:,1));
        s1_task=find(s(:,2));
        s1_subtask=find(s(:,3));
        

        if length(s1_rest)==3 && length(s1_subtask)==60
            for count=1:3
                restMat(count,1)=s1_rest(count);
                restMat(count,2)=s1_task(find(s1_task>s1_rest(count),1,'first'));
                diffTime=t(restMat(count,2))-t(s1_rest(count));
                if diffTime <= 50 || diffTime >= 70 %arbitrary boundaries, check if rest block is roughly 60 sec
                    disp(['Rest block duration is off. Check pair ' pairs{i} ' manually.']);
                    flag_hopelesscase=true;
                end
                if count<3 %take column 3 as marker for task (first subtask starting)
                    taskMat(count,1)=s1_subtask(find(s1_subtask>restMat(count,2),1,'first')); %find first subtask marker after task has officially started
                    taskMat(count,2)=s1_subtask(find(s1_subtask<s1_rest(count+1),1,'last')); %find last subtask marker before rest
                end
            end
            subtaskMat(:,1)=s1_subtask(1:30);
            subtaskMat(:,2)=s1_subtask(31:60);            
            N_done=N_done+1;
        else
            errorcount=errorcount+1;
            flag_subj2needed=true;
            disp(['subj 2 needed for pair ' pairs{i}]);
            flag_hopelesscase=true;
        end
        
        
    %load sub 2, check if different from sub 1   
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub2_' conds{c} '.mat']);
        s2_rest=find(s(:,1));
        s2_task=find(s(:,2));
        s2_subtask=find(s(:,3));
        
        if length(s1_rest)~=length(s2_rest) || length(s1_task)~=length(s2_task) || length(s1_subtask)~=length(s2_subtask) 
            if flag_subj2needed && length(s2_rest)==3 && length(s2_subtask)==60 %trying to get sensible data from subject 2 now

                for count=1:3
                    restMat(count,1)=s2_rest(count);
                    restMat(count,2)=s2_task(find(s2_task>s2_rest(count),1,'first'));
                    diffTime=t(restMat(count,2))-t(s2_rest(count));
                    if diffTime <= 50 || diffTime >= 70 %arbitrary boundaries, check if rest block is roughly 60 sec
                        disp(['Rest block duration is off. Check pair ' pairs{i} ' manually.']);
                        flag_hopelesscase=true;
                    else
                        flag_hopelesscase=false;
                    end
                    if count<3 %take column 3 as marker for task (first subtask starting)
                        taskMat(count,1)=s2_subtask(find(s2_subtask>restMat(count,2),1,'first')); %find first subtask marker after task has officially started
                        taskMat(count,2)=s2_subtask(find(s2_subtask<s2_rest(count+1),1,'last')); %find last subtask marker before rest
                    end
                end
            subtaskMat(:,1)=s2_subtask(1:30);
            subtaskMat(:,2)=s2_subtask(31:60);                 
            N_done=N_done+1;
            else
                disp('subject1 and subject2 have different number of codings.');
                flag_hopelesscase=true;
            end
        elseif any(find(s1_rest~=s2_rest,1)>0) || any(find(s1_task~=s2_task,1)>0) || any(find(s1_subtask~=s2_subtask,1)>0) %if subjects have different codings, we are having a bunch of other problems.
            disp('subject1 and subject2 have different codings');
            flag_hopelesscase=true;
        end
        
        if flag_hopelesscase
           disp(['Pair ' pairs{i} ' needs to be checked manually.']); 
        else
            %save data
            filename=[path '/Data_' conds{c} '/RPS_' pairs{i} '_marker_' conds{c} '.mat'];
            save(filename,'restMat', 'taskMat', 'subtaskMat');
        end

    end
    
