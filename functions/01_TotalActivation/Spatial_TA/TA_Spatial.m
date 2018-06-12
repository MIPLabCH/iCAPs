%% This function computes the TV regularization:
% F(x) = min ||y - x ||^2 + lambda * ||TV{x}||_1
% using the 3D extension of the FISTA:
% Beck, A., & Teboulle, M. (2009). A fast iterative shrinkage-thresholding 
% algorithm for linear inverse problems. 
% SIAM journal on imaging sciences, 2(1), 183-202.
%
% Inputs:
% y : 3D or 4D volume; if y is 4D, a 3D algorithm is applied to the first 3 
% components in parallel
% param : Structure containing relevant TA parameters; here, we require the
% fields 'NitSpat' (number of iterations for spatial regularization),
% 'LambdaSpat' (empirically tuned regularization coefficient for spatial
% regularization), 'Lip' (Lipshitz constant), 'tol' (tolerance
% threshold below which we stop iterating), 'VoxelIdx' (n_ret_vox x 3 
% coordinates of the 3D voxels used in TA) and 'Dimension' (X, Y, Z and T
% sizes)
%
% Outputs:
% - x_out : denoised output
%
% Implemented by Younes Farouj, 10.03.2016
function x_out = TA_Spatial(y,param)

    if ~isfield(param, 'NitSpat'), param.NitSpat = 400; end
    if ~isfield(param, 'LambdaSpat'), param.LambdaSpat = 2; end
    if ~isfield(param, 'Lip'), param.Lip = 12; end
    if ~isfield(param, 'tol'), param.tol = 1e-6; end

    if (length(param.Dimension) ~= 4)
        error('SIZE should have 4 dimensions.')
    end

    % param.weight is 1 if neighbouring elements are from the same tissue type,
    % and less than 1 otherwise
    if ~isfield(param, 'weight_x'), param.weight_x = ones(param.Dimension); end
    if ~isfield(param, 'weight_y'), param.weight_y =  ones(param.Dimension); end
    if ~isfield(param, 'weight_z'), param.weight_z =  ones(param.Dimension); end
    
    % This is going to be the output of the spatial regularization problem
    x_out = zeros(size(y));

    %%%% convert to 4D
    %x_vol = zeros(param.Dimension);
    y_vol = zeros(param.Dimension);
    % temp = zeros(param.Dimension);

    for i = 1:length(param.VoxelIdx(:,1));
        y_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:) = y(:,i);
    end

    % Gradient operator
    Op  =  @(x)gradient3D_full(x,param.weight_x,param.weight_y,param.weight_z);  
    
    % Its Adjoint : minus divergence <f,grad g> = -<div f,u>
    Adj_Op = @(u,v,w)-div3D_full(u,v,w,param.weight_x,param.weight_y,param.weight_z);    
    
    % compute sum TV{y}
    evaluate_norm = @(y)evaluate_3D_TV(y);   

    x_vol = MyProx(y_vol,Op,Adj_Op,evaluate_norm,param);

    for i=1:length(param.VoxelIdx(:,1));
        x_out(:,i) = x_vol(param.VoxelIdx(i,1),param.VoxelIdx(i,2),param.VoxelIdx(i,3),:);
    end

end
