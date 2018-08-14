function [iCAPs_consensus,iCAPs_nItems]=getClusterConsensus(IDX,Consensus)

% compute the average consensus of every cluster, (m(k) in Monti et al.,
% 2003)

nClus=max(IDX);
iCAPs_consensus=zeros(nClus,1);
iCAPs_nItems=zeros(nClus,1);

for iC=1:nClus
    clusID=find(IDX==iC);
    
    % ID of the sub-consensus matrix of cluster iC
    rowID=repmat(clusID,length(clusID),1);
    colID=repelem(clusID,length(clusID));
    
    subMat=zeros(size(Consensus));
    subMat(rowID,colID)=1;
    
    iCAPs_consensus(iC,1)=mean(mean(Consensus(subMat>0)));
    iCAPs_nItems(iC,1)=length(clusID);
end
