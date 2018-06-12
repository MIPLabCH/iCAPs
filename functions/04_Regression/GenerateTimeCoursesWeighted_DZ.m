%% This function computes the time courses of activity for all the assessed iCAPs.
%
% spatio-temporal regression function
% This is an alternative to obtain the time courses, the function uses the 
% information of significant innovations to costruct the design matrix.
%
% Inputs:
%     - iCAPsResults: struct containing the necessary input data
%        .iCAPs - iCAPs maps resluting from k-means clustering (n_iCAPs x n_voxels)
%        .AI -  Activity inducing signals: n_voxels x (n_subj*n_TP) matrix of 
%               concatenated subject activity-inducing data
%        .IDX - clustering information of significant innovations frames (n_significant_innov x 1)
%        .time_labels - timing information of significant innovation frames (n_significant_innov x 1)
%        .subject_labels - subject information of significant innovation frames (n_significant_innov x 1)
%        [.dist_to_centroid] - n_significant_innov X n_iCAPs matrix with 
%               distances from all significant innovation frames to all 
%               cluster centroids (iCAPs)
%      for regression with a soft assignment factor, either
%        .dist_to_centroid, or .iCAPs_folds, or .Data are required
%        [.iCAPs_folds] - information on all clustering folds, will be used to retriev distances to cluster center:
%               .sum_dist - subfield containing the total sum of distances
%                   from all frames to theirs corresponding cluster centers
%               .dist_to_centroid {} - matrix containing distances to all
%                   cluster centers for every frame
%        [.Data] - n_voxels x n_significant_innov matrix containing all 
%               significant innovation frames; required for distance to 
%               cluster computation; requires also that the field 
%               param.dystType is specified
%        [.scrub_labels] - scrubbing information of concatenated subject time courses (n_subj*n_TP x 1)
%               - will be set to all ones if not specified
%     - param: structure containing necessary parameters
%        .n_subjects - number of subjects
%        [.softClusterThres] - soft assignment factor, should be larger than
%               one (1.1 to 1.25 seem to be good values, for parameter 
%               optimization run main_getNinnovStats)
%               - if the field does not exist or is empty, hard clustering 
%                 will be used for back-projection
%        [.excludeMotionFrames] - field to select whether motion frames
%               should be excluded from the regression (if scrubbing)
%        [.n_folds] - required if distances should be taken from
%               iCAPsResults.iCAPs_folds (see above)
%
% Outputs:
%     - TC: cell array with time courses for every subject
%     - stats: structure containing model fitting statistics (RSS, BIC, AIC,...)
%
% v2.0 DZ 27.10.2017 - added compatibility with scrubbed data


function [TC,stats] = GenerateTimeCoursesWeighted(clusteringResults,param)

    %% getting variables
    iCAPs=clusteringResults.iCAPs;
    IDX=clusteringResults.IDX;
    AI=clusteringResults.AI;
    if size(AI,1)<size(AI,2)
%         warning('Inverting activity inducing signal for transient-informed regression');
        AI=AI';
    end
    
    AI_subject_labels=clusteringResults.AI_subject_labels;
    subject_labels=clusteringResults.subject_labels;
    time_labels=clusteringResults.time_labels;
    % adding option of scrubbed time courses
    if isfield(param,'excludeMotionFrames') && param.excludeMotionFrames
        scrub_labels=clusteringResults.scrub_labels;
    else
        scrub_labels=ones(size(AI,2),1); % if there was no scrubbing, include all frames in the temporal mask
    end
    
    

    
    %% constants
    nSub = param.n_subjects;
    nClus = size(iCAPs,1);
    nTP_all = size(AI,2);
    for iS=1:nSub
        nTP_sub(iS,1)=nnz(AI_subject_labels==iS);
        vols_AI_iS=AI_subject_labels==iS;
        nTP_sub_scrub(iS,1) = nnz(scrub_labels(vols_AI_iS));
    end
    nVox = size(AI,1);
    
    % get matrix with distances to cluster centroids for all frames
    if isfield(clusteringResults,'dist_to_centroid')
        dist_to_centroid=clusteringResults.dist_to_centroid;
    elseif isfield(clusteringResults,'iCAPs_folds')
        for iFold= 1:param.n_folds
            clusteringResults.iCAPs_folds.total_dist_sum(iFold,1)=sum(clusteringResults.iCAPs_folds.sum_dist{iFold});
        end
        [~,bestID]=min(clusteringResults.iCAPs_folds.total_dist_sum);
        dist_to_centroid=clusteringResults.iCAPs_folds.dist_to_centroid{bestID};
    elseif isfield(clusteringResults,'I_sig')
        dist_to_centroid=pdist2(clusteringResults.I_sig,iCAPs,param.distType);
    else
        error('For spatio-temporal regression with soft assignment factor, either the distances to the centers or the innovation frames (Data) are needed');
    end
        
    [IDX_mat,~,~]=getIDXmat(dist_to_centroid,param.softClusterThres);
    
    clear clusteringResults

    
    % compute real frame indices (without considering positive or negative
    % frames separately)
    for iS=1:nSub
        vols_iS=subject_labels==iS;
        time_labels(vols_iS&time_labels>nTP_sub(iS))=time_labels(vols_iS&time_labels>nTP_sub(iS))-nTP_sub(iS);
    end
    

%         disp('Normalizing activity inducing signals');
%         Activity_inducing_norm=zeros(size(Activity_inducing));
%         for iS=1:nSub
%             % indices of time series from subject iS
%             vols_iS=AI_subject_labels==iS;
% 
%             % divide by standard deviation
%             Activity_inducing_norm(:,vols_iS)=Activity_inducing(:,vols_iS)./...
%                 repmat(sqrt(mean(Activity_inducing(:,vols_iS).^2,2)),[1,nTP_sub]);
% 
%             % remove mean across voxels (global signal regression)
%     %         Activity_inducing_norm(:,vols_iS)=Activity_inducing_norm(:,vols_iS) - repmat(mean(Activity_inducing_norm(:,vols_iS),1),n_vox,1);
%         end
%         Activity_inducing_norm(isnan(Activity_inducing_norm))=0;
%         AI_concatenated=reshape(Activity_inducing_norm,nTP_all,[]);
%         clear Activity_inducing Activity_inducing_norm
        
    
    %% do regression for every subject
    tic
    for iS=1:nSub
        vols_AI_iS=AI_subject_labels==iS;
        % concatenate activity inducing signals
        AI_concatenated=reshape(AI(:,vols_AI_iS),nTP_sub(iS)*nVox,1); % contains concatenated AI per subject, scrubbed frame selection later
        
        tempMask_iS=scrub_labels(vols_AI_iS);
        [num_seg,startID_seg,endID_seg]=getScrubInfo(tempMask_iS);

        if num_seg>1
            fprintf('.');
        end
        
        % construct temporal mask for concatenated AI
        tempMask_iS_AI=kron(tempMask_iS,ones(nVox,1));
        
        disp(['subject: ' num2str(iS)])
        disp('Constructing design matrix ...')
        
        %number of significant innovations + start and end frame + number of disruptions due tu motion
        nInnov(iS)=nnz(subject_labels==iS)+2*nClus+(num_seg-1)*nClus; 

        S=sparse(nTP_sub_scrub(iS)*nVox,nInnov(iS));
        B=zeros(nTP_sub_scrub(iS),nInnov(iS));
        clusStartID=1;
        for iC=1:nClus
            sub_clus_time_labels{iS,iC}=sort(time_labels(subject_labels==iS&IDX_mat(:,iC)==1));
            
            % add the first frame after a motion disruption as
            % innovation
            sub_clus_time_labels{iS,iC}=[sub_clus_time_labels{iS,iC};startID_seg];
            sub_clus_time_labels{iS,iC}=unique(sub_clus_time_labels{iS,iC}); % remove doubles

            % construct new innovation time course with added
            % innovations
            newInnovTC=zeros(nTP_sub(iS),1);
            newInnovTC(sub_clus_time_labels{iS,iC})=1;

            % delete motion segments
            scrubInnovTC=newInnovTC(tempMask_iS~=0);
            if length(scrubInnovTC)~=nTP_sub_scrub(iS); error('wrong number of time points'); end

            % innovation labels to use for design matrix construction
            sub_clus_time_labels_scrub{iS,iC}=find(scrubInnovTC);
            sub_clus_time_labels_scrub{iS,iC}=[1;sub_clus_time_labels_scrub{iS,iC};length(scrubInnovTC)+1]; % add an innovation at the start and at the end
            sub_clus_time_labels_scrub{iS,iC}=unique(sub_clus_time_labels_scrub{iS,iC});

            % number of innovations per iCAP, including the first and
            % last frame and the interruptions due to excessive motion
            nInnov_Clus=length(sub_clus_time_labels_scrub{iS,iC});

            B_k=zeros(nTP_sub_scrub(iS),nInnov_Clus-1);
            for iFrame=1:nInnov_Clus-1
                firstID=sub_clus_time_labels_scrub{iS,iC}(iFrame);
                lastID=sub_clus_time_labels_scrub{iS,iC}(iFrame+1)-1;

                B_k(firstID:lastID,iFrame)=1;
            end
            S_k=sparse(kron(B_k,iCAPs(iC,:)'));
            
            S(:,clusStartID:clusStartID+nInnov_Clus-2)=S_k;
            B(:,clusStartID:clusStartID+nInnov_Clus-2)=B_k;
            
            clustID{iS,1}(clusStartID:clusStartID+nInnov_Clus-2,1)=iC;
            clusStartID=clusStartID+nInnov_Clus-1;
        end
        S(:,~sum(S,1))=[];
        nBeta=size(S,2);
        nInnov(iS)=nBeta;
        x0 = zeros(nBeta,1);
        
        disp('Solving OLS ...')
        % sparse implementation of the pseudo-inverse ((X^T X)^(-1) X^T)
        X1=S'*S; % X1=(X^T X);
        tmp=S'*AI_concatenated(tempMask_iS_AI~=0); % y1 = X^T y
        innovWeights=X1\tmp; % y2 = (X^T X)^(-1) y1
        clear X1 tmp
        
        % reconstruct time courses
        for iC=1:nClus
            betaInClusID=find(clustID{iS,1}==iC);
            innov_Clus=sub_clus_time_labels_scrub{iS,iC};
            nInnov_Clus=length(innov_Clus);
            if nInnov_Clus~=length(betaInClusID)+1
                error('there should be as many betas as piecewise constant activity sequences');
            end
            for iFrame=1:nInnov_Clus-1
                firstID=innov_Clus(iFrame);
                lastID=innov_Clus(iFrame+1)-1;
                TC_scrub{iS}(iC,firstID:lastID)=innovWeights(betaInClusID(iFrame));
            end
        end
        TC{iS}=zeros(nClus,nTP_sub(iS));
        TC{iS}(:,tempMask_iS~=0)=TC_scrub{iS};
        
        % saving mean residuals and other quality criteria (BIC, AIC,
        % ...)
        stats.RSS(iS,1)=sum((S*innovWeights-AI_concatenated(tempMask_iS_AI~=0)).^2); % residual sum of squares
        stats.n(iS,1)=size(S,1); % number of observations (voxels*frames)
        stats.k(iS,1)=size(S,2); % number of regressors (amplitude segments to estimate)
        stats.bic(iS,1)=stats.n(iS)*log(stats.RSS(iS)/stats.n(iS)) + stats.k(iS)*log(stats.n(iS));
        stats.aic(iS,1)=stats.n(iS)*log(stats.RSS(iS)/stats.n(iS)) + stats.k(iS)*2;
    end
    toc
end