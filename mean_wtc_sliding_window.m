%sliding window with increments of 1 sample
%
%-wtc_signal needs to be signal in frequency range of interest, and of
%length of the trial.
%-realtime_cutoff is calculated cutoff of cone of influence as if signal is
%measured in realtime (put '0' if no pseudo-real time) in samples.
%-windowsize is desired size of window in samples.


function [mean_wtc_window_vec]=mean_wtc_sliding_window(wtc_signal, realtime_cutoff, windowsize)

            mean_wtc_window_vec=NaN(length(wtc_signal)-windowsize,1);
            
            for w=1:length(wtc_signal)-windowsize
                
                teststat=wtc_signal(w:w+windowsize);
                
                for l=1:length(teststat)
                    if l<realtime_cutoff || l>length(teststat)-realtime_cutoff
                        teststat(l)=NaN;
                    end
                end
                mean_wtc_window_vec(w)=mean(teststat,'omitnan'); 
            end
                
                
            


%             for v=windowlength:T*triallength%length(xvec)%-windowlength %T*triallength=total number of samples per block
%             %first block
%             start=subtaskMat(1,1);   
%          %   disp(num2str(start));
%                assert(start+v<=restMat(2,1), ['sliding window running out of first task block, please check. start+v=' num2str(start+v-1) ', rest=' num2str(restMat(2,1))]);
%                     teststat=soi(start+v-windowlength:start+v-1);
%             
%             for l=1:length(teststat)
%                 if l<cutoff || l>length(teststat)-cutoff
%                     teststat(l)=NaN;
%                 end
%             end
%             mean_wtc_window_vec(v)=mean(teststat,'omitnan'); 
%             
%             end
            
end