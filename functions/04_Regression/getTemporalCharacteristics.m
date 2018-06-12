%% This function computes the temporal characteristics (duration, occurences, etc.)
% given a set of time courses
%
% Inputs:
% - TC: cell array with iCAPs time courses for every subject
% - clusteringResults: structure containing clustering results, fields:
%       .subject_labels
%       .time_labels
%       [.scrub_labels]
% - param: structure with input parameters
%       .TR: TR
%
% Outputs:
% - out: structure containing all computed characteristics:
%       .durations_total_n_subs: total number of frames an iCAP is active 
%                   in every subject
%       .durations_total_s_subs: total time (in seconds) an iCAP is active 
%                   in every subject
%
%       .durations_avg(_n/_s): average duration of an iCAPs activity

function [tempChar] = getTemporalCharacteristics(TC,clusteringResults,param)
    
    nSub=length(TC);
    for iS=1:nSub
        nTP_sub(iS)=nnz(clusteringResults.AI_subject_labels==iS);
    end
    subject_labels=clusteringResults.subject_labels;
    IDX=clusteringResults.IDX;
    
    % adding option of scrubbed time courses
    if isfield(param,'excludeMotionFrames') && param.excludeMotionFrames
        scrub_labels=iCAPsResults.scrub_labels;
    else
        scrub_labels=ones(nTP_sub*nSub,1);
    end

    
        tempChar = computeTemporalCharacteristics(iCAPsResults.Time_Courses,1,param,scrub_labels,subject_labels,IDX);
%         plotTCandInnov(iCAPsResults,param);


end


%%
function plotTCandInnov(iCAPsResults,param)
    %% constants
    nSub=length(iCAPsResults.Time_Courses);
    nClus=size(iCAPsResults.Time_Courses{1},1);
    nTP_sub=size(iCAPsResults.Time_Courses{1},2);
    
    % adding option of scrubbed time courses
    if isfield(param,'excludeMotionFrames') && param.excludeMotionFrames
        scrub_labels=iCAPsResults.scrub_labels;
    else
        scrub_labels=ones(nTP_sub*nSub,1);
    end

    for iS=1:nSub
        vols_iS=(iS-1)*nTP_sub+1:iS*nTP_sub;
        nTP_sub_scrub(iS,1) = nnz(scrub_labels(vols_iS));
        TemporalMask{iS}=scrub_labels(vols_iS);
    end

    
    %% get data and iCAPs output directories
    if isfield(param,'data_title') % in case of new data saving structure (clustering in subfolders)
        % there is a main folder according to data+thresholding
        % and a separate one for clustering
        outDir_main=fullfile(param.PathData,'iCAPs_results',[param.data_title,'_',param.thresh_title]);
        outDir_iCAPs=fullfile(param.PathData,'iCAPs_results',[param.data_title,'_',param.thresh_title],param.iCAPs_title);
    else
        % before, there was a separate folder for
        % thresholding/data/clustering parameters
        outDir_main=fullfile(param.PathData,'iCAPs_results',[param.thresh_title,'_',param.iCAPs_title]);
        outDir_iCAPs=fullfile(param.PathData,'iCAPs_results',[param.thresh_title,'_',param.iCAPs_title]);
    end

    
    %% plotting time courses and innovation indices
    if ~exist(fullfile(outDir_iCAPs,'Time_Courses_and_Innovations'),'dir'); mkdir(fullfile(outDir_iCAPs,'Time_Courses_and_Innovations'));end;
    if ~exist(fullfile(outDir_iCAPs,'Time_Courses'),'dir'); mkdir(fullfile(outDir_iCAPs,'Time_Courses'));end;
    
    iCAPsResults.time_labels(iCAPsResults.time_labels>195)=-(iCAPsResults.time_labels(iCAPsResults.time_labels>195)-195);

    for iS = 1:nSub
        % normalize time courses
        iCAPsResults.Time_Courses_plot{iS,1}=reshape(zscore(iCAPsResults.Time_Courses{iS}(:)),size(iCAPsResults.Time_Courses{iS},1),size(iCAPsResults.Time_Courses{iS},2));
    end
    
    
    colors = cbrewer('div', 'Spectral', 10);
    cmap=cbrewer('div', 'RdYlBu',100);
    cmap=cmap(100:-1:1,:);
    %     colors = [colors;colors([11:-1:1],:)];
    set(groot,'defaultAxesColorOrder',colors);
    set(groot,'defaultAxesFontSize',15);

    clims=[-3,3];
    
    for iS=1:10:nSub
        figure('position',[440   378   515   420]);
        imagesc(iCAPsResults.Time_Courses_plot{iS},clims);
        title(['subject ' num2str(iS)])
        colormap(cmap)
        c=colorbar;
        c.Label.String='amplitude';
        xlabel('time [frames]');
        % ylabel('iCAP');
        print(fullfile(outDir_iCAPs,'Time_Courses',['SUB' num2str(iS)]),'-dpng','-painters');

        
        for iC=1:nClus
            fig=figure('Position',[440   573   560   167]);
            hold on;
            grid on
            stemAmp=max(abs(iCAPsResults.Time_Courses_plot{iS,1}(iC,:)))/2;
            stemAmp=2;
            patch('Vertices',[0 -1; nTP_sub -1; nTP_sub 1; 0 1],'Faces',[1 2 3 4],'facecolor',0.5*[1 1 1],'edgeColor',0.5*[1 1 1],'FaceAlpha',0.3,'EdgeAlpha',0.3);
            pl(1)=stem(abs(iCAPsResults.time_labels(iCAPsResults.subject_labels==iS&iCAPsResults.IDX==iC)),...
                stemAmp*sign(iCAPsResults.time_labels(iCAPsResults.subject_labels==iS&iCAPsResults.IDX==iC)),'k')
            plot([0 200],[0 0],'k-')
            pl(2)=plot(iCAPsResults.Time_Courses_plot{iS}(iC,:),'linewidth',1.5,'color',colors(3,:));
    %         title(['subject ' num2str(iS)])
            xlabel('time [frames]')
    %         ylabel('amplitude')
    %             set(gca,'ylim',2.5*[-stemAmp stemAmp]);
            set(gca,'ytick',[-20:2:20],'ylim',[min(-4,1.1*min(iCAPsResults.Time_Courses_plot{iS,1}(iC,:))) max(4,1.1*max(iCAPsResults.Time_Courses_plot{iS,1}(iC,:)))]);
            if ~exist(fullfile(outDir_iCAPs,'Time_Courses_and_Innovations'),'dir'); mkdir(fullfile(outDir,'Time_Courses_and_Innovations'));end;
            outFileName=fullfile(outDir_iCAPs,'Time_Courses_and_Innovations',['SUB' num2str(iS) ' - iCAP' num2str(iC)]);
    %             figure(fig);
    %         legend(pl,{'',''},'location','eastoutside');
            set(gca,'xlim',[1 195])
            print(outFileName,'-depsc2','-painters');
            close gcf
            
            
            
%             stemAmp=max(abs(iCAPsResults.Time_Courses_plot{iS}(iC,:)))/2;
%             fig=figure('position',[440   662   560   136]);
%             imagesc([1 nTP_sub]+1,[-5 5],TemporalMask{iS}','alphadata',0.1,[0,1]);colormap('gray');
%             hold on;
%             stem(abs(iCAPsResults.time_labels(iCAPsResults.subject_labels==iS&iCAPsResults.IDX==iC)),...
%                 stemAmp*sign(iCAPsResults.time_labels(iCAPsResults.subject_labels==iS&iCAPsResults.IDX==iC)))
%             stairs(iCAPsResults.Time_Courses_plot{iS}(iC,:),'linewidth',1.5);
%             title(['SUB' num2str(iS) ' - iCAP' num2str(iC)])
%             xlabel('frame')
% %             ylabel('amplitude')
%             set(gca,'ylim',2.5*[-stemAmp stemAmp]);
%             outFileName=fullfile(outDir_iCAPs,'Time_Courses_and_Innovations',['SUB' num2str(iS) ' - iCAP' num2str(iC)]);
%             figure(fig);
%             print(outFileName,'-depsc2','-painters');
%             close gcf
        end
    end
end
