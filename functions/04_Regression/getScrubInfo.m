function [num_seg,startID_seg,endID_seg] = getScrubInfo(TM)


% compute the number of segments that are together
diff_TM=[0;diff(TM)]; % first TP is zero to fit to same length as TM

if length(find(diff_TM==-1))==length(find(diff_TM==1))
    if ~TM(1)&&~TM(end) % if moved at beginning AND at the end
        num_seg=length(find(diff_TM==-1));
        startID_seg=find(diff_TM==1);
        endID_seg=find(diff_TM==-1)-1;
    else % if not moved neither at beginning nor at end
        num_seg=length(find(diff_TM==-1))+1;
        startID_seg=[1;find(diff_TM==1)]; 
        endID_seg=[find(diff_TM==-1)-1;length(TM)];
    end
else % if the number of activations and deactivations is not the same, one block is at the end or the beginning ot the signal
    num_seg=max(length(find(diff_TM==-1)),length(find(diff_TM==1)));
    
    if ~TM(1) % moved at start
        startID_seg=find(diff_TM==1);
        endID_seg=[find(diff_TM==-1)-1;length(TM)];
    end
    
    if ~TM(end) % moved at end
        startID_seg=[1;find(diff_TM==1)];
        endID_seg=find(diff_TM==-1)-1;
    end
end