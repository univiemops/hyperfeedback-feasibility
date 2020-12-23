%get markers from files 's' for 'rest' and 'task' (overall task)

%KK 14.12.20
%assume three rest blocks and two task blocks per file

%assume furthermore: (valid at least for FP)
%s(:,1) is rest
%s(:,2) is task
%s(:,3) is also task, more specific task. Called 'subtask' furtheronwards.
%s always demarks starting point of incidents.

% Seems like for PS it holds: s(:,8) is rest, s(:,9) is task, s(:,10) is
% some subtask...

%for rest take beginning of rest to task marker
%for task take first and last subtask marker


%loop over all subjects and conditions
%load in 's'
%calculate diffs
%write into matrices to safe later on

%check no-matches manually for safety!



path=['/Users/kathrin/data/RPS/Kathrin_collab/hmrData'];

pairs={'01','02','03','04','05','06','08','09','10','11','13','14','15','17','18','19','20','21','22','23','24','26','28','29','30','31','32'};
conds={'FP','PS', 'PD','C'};
flag_nRestgrtrnTask=false;

nomatch=zeros(length(pairs),length(conds))*NaN;
%Nch=16; %Number of channels



for c=1%:length(conds)
    
    pairmat=zeros(length(pairs),1)*NaN;
    
    errorcount=0;
    errorcount_correctRandT=0;
    errorcount_correctRnotT=0;
    N_done=0;
    
    for i=1:length(pairs) %i=15: pair 18 is a faultless pair
        flag_nRestgrtrnTask=false;
        flag_subj2needed=false;
        flag_hopelesscase=false;
        
        restMat=zeros(3,2)*NaN;
        taskMat=zeros(2,2)*NaN;
        
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub1_' conds{c} '.mat']);
        %get all the markers for each condition
        s1_rest=find(s(:,1));
        s1_task=find(s(:,2));
        s1_subtask=find(s(:,3));
        
        deltaRestTask=zeros(length(s1_rest),1);
        %   deltaRestSubtask=zeros(length(s1_rest),1);
        
        
%         %case correct rest and task codings
%         if length(s1_rest)==3 && length(s1_task)==3 %assuming three times rest block, ending with task marker
%             disp('hooray, assuming standard case');
%             for count=1:3
%                 if s1_rest(count)< s1_task(count)
%                     restMat(count,1)=s1_rest(count);
%                     restMat(count,2)=s1_task(count);
%                     diffTime=t(s1_task(count))-t(s1_rest(count));
%                     
%                     if diffTime <= 50 || diffTime >= 70 %arbitrary boundaries, check if rest block is roughly 60 sec 
%                        disp(['Rest block duration is off.']);
%                        flag_hopelesscase=true;
%                     end
%                     if count<3 %take column 3 as marker for task (first subtask starting)
%                         taskMat(count,1)=s1_subtask(find(s1_subtask>s1_task(count),1,'first')); %find first subtask marker after task has officially started
%                         taskMat(count,2)=s1_subtask(find(s1_subtask<s1_rest(count+1),1,'last')); %find last subtask marker before rest
%                     end
%                     
%                 else %I have to think
%                     disp(['markers in subject 1 of pair ' pairs{i} 'are in correct number but wrong order. Trying out if subject 2 is ok.']);
%                     flag_subj2needed=true;
%                     continue
%                 end
%             end
%             N_done=N_done+1;
%         else
%             errorcount_correctRandT=errorcount_correctRandT+1;
%             flag_subj2needed=true;
%             %      disp(['Houston, we are having a problem. Pair ' pairs{i} ' subj 1 does not have the correct amount of markers.']);
%         end
        
        %case correct rest codings, trying to work around wrong amount of task
        %codings
        if length(s1_rest)==3 %&& length(s1_task)~=3
        %    disp('other case');
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
            N_done=N_done+1;
        else
            % errorcount_correctRnotT=errorcount_correctRnotT+1;
            errorcount=errorcount+1;
            flag_subj2needed=true;
            disp(['subj 2 needed for pair ' pairs{i}]);
            flag_hopelesscase=true;
        end
        
        
        
        
        
        % %case rescuing stuff with wrong codings
        %         if length(s1_task)<length(s1_rest) && flag_subj2needed==false
        %             disp(['more rest than task markers for pair ' pairs{i}]);
        %             flag_nRestgrtrnTask=true;
        %         end
        %         %  if length(s1_subtask)<length(s1_rest), error('more rest than subtask markers'); end %noch nicht das End vom Lied: du kannst mehr Marker haben, die trotzdem vor dem letzten rest-marker liegen....
        %         if s1_task(end)<s1_rest(end)
        %             disp(['rest block has no end marker for subject 1 of pair ' pairs{i} '. Skipping this pair.']);
        %             continue
        %         end
        %
        %
        %         %find all differences to major task markers/blocks of rest
        %         index=1;
        %         while index<=length(s1_rest) && flag_subj2needed==false
        %             taskindex=find(s1_task>s1_rest(index),1,'first');
        %             if flag_nRestgrtrnTask == false
        %                 deltaRestTask(index)=(s1_task(taskindex)-s1_rest(index))/fs; %fs sampling rate
        %             elseif flag_nRestgrtrnTask == true && index<length(s1_rest)
        %                 if s1_rest(index+1)< s1_task(taskindex)
        %                     deltaRestTask(index)=(s1_task(taskindex)-s1_rest(index+1))/fs;
        %                     index=index+1;
        %                 else
        %                     deltaRestTask(index)=(s1_task(taskindex)-s1_rest(index))/fs;
        %                 end
        %             end
        %             index=index+1;
        %         end
        %
        % %         for index=1:length(s1_rest) %find all differences to major task markers/blocks of rest
        % %         taskindex=find(s1_subtask>s1_rest(index),1,'first');
        % %         deltaRestSubtask(index)=(s1_subtask(taskindex)-s1_rest(index))/fs; %fs sampling rate
        % %         end
        %
        load([path '/Data_' conds{c} '/RPS_' pairs{i} '_sub2_' conds{c} '.mat']);
        s2_rest=find(s(:,1));
        s2_task=find(s(:,2));
        s2_subtask=find(s(:,3));
        
        if length(s1_rest)~=length(s2_rest) || length(s1_task)~=length(s2_task) || length(s1_subtask)~=length(s2_subtask)
            if flag_subj2needed && length(s2_rest)==3 %trying to get sensible data from subject 2 now
                % if length(s1_rest)==3 %&& length(s1_task)~=3 %I think this could be general case, need to check
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
            %safe data
            filename=[path '/Data_' conds{c} '/RPS_' pairs{i} '_marker_' conds{c} '.mat'];
            save(filename,'restMat', 'taskMat');
        end
        %
        %         if flag_subj2needed %if stuff does not work with subj1 try if at least subj2 has useful data
        %             %dummy
        %         end
        %
        %
        %         pairmat(i,1)=mean(deltaRestTask);
    end
    
end