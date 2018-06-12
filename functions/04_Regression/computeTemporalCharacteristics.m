% computes temporal characteristics of iCAPs
% occurrences of only one frame will not be counted


% Input: TC - nSubjects x 1 cell object
%             cells: nClus x nTP time courses per iCAP
%        clusteringResults - strcut with clustering results
%           .AI_subject_labels
%           .subject_labels
%           .IDX
%           [.scrub_labels] - only needed if param.excludeMotionFrames=1
%        param - necessary field: 
%             .activityThres - threshold at which normalized time courses 
%               will be considered "active", default = 1
%        
%
% Output: tempChar - structure containing fields with temporal
%                    characteristics (nClus x nSubs)
%       * Characteristics of significant innovations
%           .innov_counts_total - number of all significant innovations per
%               subject
%           .innov_counts_total_perc - percentage of significant 
%               innovations in all included time points
%           .innov_counts - nClus x nSub matrix, number of significant
%               innovations per iCAP per subject
%           .innov_counts_percOfInnov - nClus x nSub matrix, percentage of 
%               innovations in this iCAP in all significant innovations
% 
%       * Thresholded time courses and overall characteristics
%         (nSubjects x 1 cell objects)
%           .TC_norm_thes - normalized and thresholded time courses, with
%               occurences of only one frame removed from the TC
%           .TC_active - activity information, 1 if positive activity, -1
%               if negative, 0 if none
%           .coactiveiCAPs_total - 1 x nTP_subject with number of active 
%               iCAPs per TP
%           .coactivation - nClus x nClus x nTP_subject with coactivation
%               between every two iCAPs 
% 
%       * Temporal characteristics of activity blocks in time courses
%         (nClus x nSub matrices)
%           .occurences - number of activity blocks in thresholded time 
%               courses for each cluster in each subject
%           .occurences_pos - number of positive activity blocks in 
%               thresholded time courses
%           .occurences_neg - number of negative activity blocks in 
%               thresholded time courses
% 
%           .durations_total_counts - number of active time points per 
%               cluster per subject
%           .durations_total_pos_counts - number of positively active time 
%               points per cluster per subject
%           .durations_total_neg_counts - number of negatively active time 
%               points per cluster per subject
% 
%           .durations_total_perc - total duration of iCAP in percentage of
%               the whole scan duration
%           .durations_total_pos_perc - total positive duration of iCAP in 
%               percentage of the whole scan duration
%           .durations_total_neg_perc - total negative duration of iCAP in 
%               percentage of the whole scan duration
% 
%           .durations_avg_counts - average duration of activity blocks 
%               (number of time points)
%           .durations_avg_pos_counts - average duration of positive 
%               activity blocks (number of time points)
%           .durations_avg_pos_counts - average duration of negative 
%               activity blocks (number of time points)
% 
%       * Co-activation characteristics of iCAPs time courses
%         (nClus x nClus x nSub)
%           .coactivation_counts - coactivation duration (number of time
%               points)
%           .coactivation_normBoth - coactivation duration (percentage of
%               total duration of both iCAPs)
% 
%           .coactivation_sameSign - same-signed coactivation (positive
%               coupling) duration (number of time points)
%           .coactivation_diffSign - differently-signed coactivation
%               (negative coupling) duration (number of time points)
%           .coactivation_sameSign_normBoth - same-signed coactivation
%               duration (percentage of tota duration of both iCAPs)
%           .coactivation_diffSign_normBoth - differently-signed
%               coactivation duration (percentage of tota duration of both 
%               iCAPs)


function tempChar = computeTemporalCharacteristics(TC,clusteringResults,param)

thres=param.activityThres;
AI_subject_labels=clusteringResults.AI_subject_labels;
subject_labels=clusteringResults.subject_labels;
IDX=clusteringResults.IDX;

% constants
nSub=length(TC);
nClus=size(TC{1},1);
for iS=1:nSub
    nTP_sub(iS)=size(TC{iS},2);
end
nTP_all=length(AI_subject_labels);

if ~isfield(clusteringResults,'scrub_labels') || isempty(clusteringResults.scrub_labels) ...
   || ~isfield(param,'excludeMotionFrames') || ~param.excludeMotionFrames
    scrub_labels=ones(nTP_all,1);
else
    scrub_labels=clusteringResults.scrub_labels;
end


% loop over all subjects
for iS = 1:nSub
    % scrubbed time points
    vols_iS=AI_subject_labels==iS;
    tempChar.scrub_labels{iS}=scrub_labels(vols_iS)';
    tempChar.nTP_sub(iS,1)=nnz(scrub_labels(vols_iS));
    
    % innovation counts per subject
    tempChar.innov_counts_total(1,iS)=nnz(subject_labels==iS); % number of all significant innovations per subject
    tempChar.innov_counts_total_perc(1,iS)=tempChar.innov_counts_total(1,iS)/tempChar.nTP_sub(iS)*100; % percentage of significant innovations in all included time points
    
    % normalize time courses
    TC_norm{iS}=reshape(zscore(TC{iS}(:)),size(TC{iS},1),size(TC{iS},2));
    TC_norm_thes=TC_norm;
    TC_norm_thes{iS}(abs(TC_norm{iS})<thres)=0;
    
    % remove occurences of only one frame
    for iC=1:size(TC_norm_thes{iS},1) % I am computing thresholded time courses for all iCAPs first because I will use them for the co-activation computation
        activeComp=bwconncomp(TC_norm_thes{iS}(iC,:));
        
        % splitting components with sign changes
        activeComp=splitPosNegComps(activeComp,TC_norm_thes{iS}(iC,:));
        
        % deleting occurrences of only one frame
        TC_norm_thes{iS}(iC,:)=deleteOneFrameOcc(activeComp,TC_norm_thes{iS}(iC,:));
    end
    
    % number of co-active iCAPs per time point
    tempChar.TC_norm_thes{iS}=TC_norm_thes{iS};
    tempChar.TC_active{iS}=TC_norm_thes{iS};
    tempChar.TC_active{iS}(tempChar.TC_active{iS}>0)=1;
    tempChar.TC_active{iS}(tempChar.TC_active{iS}<0)=-1;
    
    tempChar.coactiveiCAPs_total{iS}=sum(tempChar.TC_active{iS}~=0);
    tempChar.coactiveiCAPs_total{iS}(tempChar.scrub_labels{iS}==0)=nan; % set not included time points to nan
    
    
    
    % compute iCAPs-wise characteristics
    for iC=1:size(TC_norm_thes{iS},1)
        % compute the active components again
        activeComp=bwconncomp(TC_norm_thes{iS}(iC,:));
        activeComp=splitPosNegComps(activeComp,TC_norm_thes{iS}(iC,:));
        
        % get signs of components
        activeComp=getCompSign(activeComp,TC_norm_thes{iS}(iC,:));

        for iA=1:activeComp.NumObjects
            if length(activeComp.PixelIdxList{iA})<2
                error('Something went wrong while suppressing occurrences of only one frame');
            end
        end
        
        % compute occurences and average durations
        tempChar.occurences(iC,iS)=activeComp.NumObjects;
        tempChar.occurences_pos(iC,iS)=nnz(activeComp.compSign>0);
        tempChar.occurences_neg(iC,iS)=nnz(activeComp.compSign<0);
        
        tempChar.durations_total_counts(iC,iS)=nnz(TC_norm_thes{iS}(iC,:));
        tempChar.durations_total_pos_counts(iC,iS)=nnz(TC_norm_thes{iS}(iC,:)>0);
        tempChar.durations_total_neg_counts(iC,iS)=nnz(TC_norm_thes{iS}(iC,:)<0);
        
        tempChar.durations_avg_counts(iC,iS)=tempChar.durations_total_counts(iC,iS)/activeComp.NumObjects;
        tempChar.durations_avg_pos_counts(iC,iS)=tempChar.durations_total_pos_counts(iC,iS)/tempChar.occurences_pos(iC,iS);
        tempChar.durations_avg_neg_counts(iC,iS)=tempChar.durations_total_neg_counts(iC,iS)/tempChar.occurences_neg(iC,iS);
        
        % innovation counts
        tempChar.innov_counts(iC,iS)=nnz(IDX==iC&subject_labels==iS);
        tempChar.innov_counts_percOfInnov(iC,iS)=tempChar.innov_counts(iC,iS)/tempChar.innov_counts_total(1,iS)*100; % percentage of innovations in this iCAP in all significant innovations        
        
%         % innovation counts in time courses
%         tempChar.innov_counts_TC(iC,iS)=nnz(diff(TC_norm{iS}(iC,:))~=0)-nnz(diff(tempChar.scrub_labels{iS}~=0)); % significant innovations in time courses, except those due to motion scrubbing
%         tempChar.innov_counts_TC_perc(iC,iS)=tempChar.innov_counts_TC(iC,iS)/tempChar.nTP_sub(iS)*100; % percentage of innovations in this iCAP in included time points per subject
%         
%         % number of co-active iCAPs with given iCAP
%         tempChar.coactiveiCAPs(iC,iS,:)=tempChar.coactiveiCAPs_total(iS,:);
%         tempChar.coactiveiCAPs(iC,iS,tempChar.TC_active{iS}(iC,:)==0)=nan; % set to nan all time points where the given iCAP was not active
%         tempChar.coactiveiCAPs_pos(iC,iS,:)=tempChar.coactiveiCAPs(iC,iS,:);
%         tempChar.coactiveiCAPs_pos(iC,iS,tempChar.TC_active{iS}(iC,:)<=0)=nan; % set to nan all time points where the given iCAP was not positive
%         tempChar.coactiveiCAPs_neg(iC,iS,:)=tempChar.coactiveiCAPs(iC,iS,:);
%         tempChar.coactiveiCAPs_neg(iC,iS,tempChar.TC_active{iS}(iC,:)>=0)=nan; % set to nan all time points where the given iCAP was not negative
%         
        
        % compute signed co-activations with all other iCAPs
        for iC2=1:nClus
            % time points of co-activation of iCAP iC and iCAP iC2
            tempChar.coactivation{iS}(iC,iC2,:)=TC_norm_thes{iS}(iC,:) & ...
                TC_norm_thes{iS}(iC2,:);
            
            % percentage of co-activation with iCAP iC2, with respect to
            % total activation of iCAP iC
            tempChar.coactivation_counts(iC,iC2,iS)=nnz(tempChar.coactivation{iS}(iC,iC2,:));
            tempChar.coactivation_normBoth(iC,iC2,iS)=nnz(tempChar.coactivation{iS}(iC,iC2,:))/...
                nnz(TC_norm_thes{iS}(iC,:)~=0 | TC_norm_thes{iS}(iC2,:)~=0);
            
            % signed co-activation
            tempChar.coactivation_posPos{iS}(iC,iC2,:)=TC_norm_thes{iS}(iC,:)>0 & ...
                TC_norm_thes{iS}(iC2,:)>0;
            tempChar.coactivation_posNeg{iS}(iC,iC2,:)=TC_norm_thes{iS}(iC,:)>0 & ...
                TC_norm_thes{iS}(iC2,:)<0;
            tempChar.coactivation_negPos{iS}(iC,iC2,:)=TC_norm_thes{iS}(iC,:)<0 & ...
                TC_norm_thes{iS}(iC2,:)>0;
            tempChar.coactivation_negNeg{iS}(iC,iC2,:)=TC_norm_thes{iS}(iC,:)<0 & ...
                TC_norm_thes{iS}(iC2,:)<0;
            
            tempChar.coactivation_sameSign(iC,iC2,iS)=(nnz(tempChar.coactivation_posPos{iS}(iC,iC2,:))+...
                nnz(tempChar.coactivation_negNeg{iS}(iC,iC2,:)));
            tempChar.coactivation_diffSign(iC,iC2,iS)=(nnz(tempChar.coactivation_posNeg{iS}(iC,iC2,:))+...
                nnz(tempChar.coactivation_negPos{iS}(iC,iC2,:)));
            
            % percentage of signed co-activation with iCAP iC2, with
            % respect to total positive or negative activation of both
            % iCAPs
            tempChar.coactivation_sameSign_normBoth(iC,iC2,iS)=(nnz(tempChar.coactivation_posPos{iS}(iC,iC2,:))+...
                nnz(tempChar.coactivation_negNeg{iS}(iC,iC2,:)))/...
                nnz(TC_norm_thes{iS}(iC,:)~=0 | TC_norm_thes{iS}(iC2,:)~=0);
            tempChar.coactivation_diffSign_normBoth(iC,iC2,iS)=(nnz(tempChar.coactivation_posNeg{iS}(iC,iC2,:))+...
                nnz(tempChar.coactivation_negPos{iS}(iC,iC2,:)))/...
                nnz(TC_norm_thes{iS}(iC,:)~=0 | TC_norm_thes{iS}(iC2,:)~=0);
            
            
            if iC==iC2;
                tempChar.coactivation_counts(iC,iC2,iS)=nan;
                tempChar.coactivation_normBoth(iC,iC2,iS)=nan;
                
                tempChar.coactivation_sameSign(iC,iC2,iS)=nan;
                tempChar.coactivation_diffSign(iC,iC2,iS)=nan;
                
                tempChar.coactivation_sameSign_normBoth(iC,iC2,iS)=nan;
                tempChar.coactivation_diffSign_normBoth(iC,iC2,iS)=nan;
            end
            
        end
        
    end
    
    % remove nans
    tempChar.durations_avg_counts(isnan(tempChar.durations_avg_counts))=0; 
    tempChar.durations_avg_pos_counts(isnan(tempChar.durations_avg_pos_counts))=0; 
    tempChar.durations_avg_neg_counts(isnan(tempChar.durations_avg_neg_counts))=0;
    
    % compute measures in percentage
    tempChar.durations_total_perc(:,iS)=tempChar.durations_total_counts(:,iS)/tempChar.nTP_sub(iS)*100;
    tempChar.durations_total_pos_perc(:,iS)=tempChar.durations_total_pos_counts(:,iS)/tempChar.nTP_sub(iS)*100;
    tempChar.durations_total_neg_perc(:,iS)=tempChar.durations_total_neg_counts(:,iS)/tempChar.nTP_sub(iS)*100;
end


% % compute measures in seconds
% tempChar.durations_total_sec=tempChar.durations_total_counts.*param.TR;
% tempChar.durations_avg_sec=tempChar.durations_avg_counts.*param.TR;
% tempChar.durations_total_pos_sec=tempChar.durations_total_pos_counts.*param.TR;
% tempChar.durations_avg_pos_sec=tempChar.durations_avg_pos_counts.*param.TR;
% tempChar.durations_total_neg_sec=tempChar.durations_total_neg_counts.*param.TR;
% tempChar.durations_avg_neg_sec=tempChar.durations_avg_neg_counts.*param.TR;
    
end




function activeComp=splitPosNegComps(activeComp,TC)
    for iA=1:activeComp.NumObjects
        % split connected components that include a sign change
        sign_comp=sign(TC(activeComp.PixelIdxList{iA}));
        sign_diff=find(diff(sign_comp));
        if ~isempty(sign_diff)
            activeComp.NumObjects=activeComp.NumObjects+length(sign_diff);
            sign_diff=[sign_diff,length(activeComp.PixelIdxList{iA})];
            % add additional components
            for iN=1:length(sign_diff)-1
                activeComp.PixelIdxList{end+1}=activeComp.PixelIdxList{iA}(sign_diff(iN)+1:sign_diff(iN+1));
            end
            % update existing component (only keep the first connected
            % component)
            activeComp.PixelIdxList{iA}=activeComp.PixelIdxList{iA}(1:sign_diff(1));
        end
    end
end


function TC=deleteOneFrameOcc(activeComp,TC)
    % deleting occurrences of only one frame
    for iA=1:activeComp.NumObjects
        if length(activeComp.PixelIdxList{iA})<2
            TC(activeComp.PixelIdxList{iA})=0;
        end
    end
end


function activeComp=getCompSign(activeComp,TC)
    activeComp.compSign=[];
    for iA=1:activeComp.NumObjects
        % split connected components that include a sign change
        sign_comp=sign(TC(activeComp.PixelIdxList{iA}));
        sign_diff=find(diff(sign_comp));
        if ~isempty(sign_diff)
            activeComp=splitPosNegComps(activeComp,TC);
            activeComp=getCompSign(activeComp,TC);
            continue;
        else
            activeComp.compSign(iA)=sign_comp(1);
        end
    end
end
