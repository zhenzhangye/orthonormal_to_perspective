function [normals, dz, N_unnormalized, zx, zy, xx, yy, Dx, Dy, J_n_un, J_dz] = Depth2Normals_o2p(z, K, mask, grad_option, xx, yy, Dx, Dy)
%depth2Normals is a function to calculate the normal map from a depth map
%based on perspective projection
%INPUT:
%       z     is the depth map of size mxn (if mask is given) or of size mnx1
%             (if no mask is given).
%       K     is the 3x3 matrix corresponding to the intrinsic parameters
%       mask  is an mxn binary mask (only necessary, if xx,yy,Dx,Dy) are not
%             known so far. If xx,yy,Dx,Dy are known, mask can be left empty
%             i.e., mask = [].
%       xx    are the pixels in x-direction wrt the principal point K(1,3)
%       yy    are the pixels in y-direction wrt the principal point K(2,3)
%       Dx    is either a sparse matrix representing the gradient in
%             x-direction or a non-sparse vector representing the gradient of z
%             in x-direction.
%       Dy    is either a sparse matrix representing the gradient in
%             y-direction or a non-sparse vector representing the gradient of z
%             in y-direction.
%OUTPUT:
%       normals are the normals of size mnx3
%       dz      is the unnormalized norm of the normals before
%               normalization. It is strongly related to the minimal surface
%               element.
%       zx      is the derivative of z in x-direction
%       zy      is the derivative of z in y-direction
%       xx      are the pixels in x-direction wrt the principal point K(1,3)
%       yy      are the pixels in y-direction wrt the principal point K(2,3)
%       Dx      is either a sparse matrix representing the gradient in
%               x-direction.
%       Dy      is either a sparse matrix representing the gradient in
%               y-direction.
%
%OPTIONAL OUTPUT:
%       J       is the Jacobian matrix J(zx,zy,z)
%       dz_p    represents K(1,1) * n_x - xx .* n_z
%       dz_q    represents K(2,2) * n_y - yy .* n_z
%       dz_z    represents -n_z
%
%EXAMPLE: You can call this function in three different ways. If you call this function
%for the first time you should call it as described in 1., after that you
%should call it as described in 2. or 3.
%         1. z is some mxn matrix representing depth values of type double
%           [N, dz, zx, zy, xx, yy, Dx, Dy] = depth2Normals(z, K, mask);
%
%         2. z is some mnx1 vector representing depth values of type double
%            xx and yy are vectors, same as the output from 1.
%            Dx, Dy are sparse matrices, same as the output from depth2Normals1.
%            [N, dz, zx, zy] = depth2Normals(z, K, [], xx, yy, Dx, Dy);
%         3. z is some mnx1 vector representing depth values of type double
%            xx and yy are vectors, same as the output from 1.
%            zx, zy are vectors vectors describing the derivative of z wrt {x,y}
%            [N, dz] = depth2Normals(z, K, [], xx, yy, zx, zy);
% Copyright by
% Author: Bjoern Haefner
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if isvector(K) % K is 2x1 or 1x2
    if length(K)==2
        cam_model = 'orthographic';
    else
        error('Error in depth2Normals(): K is a vector of length %d, but it must be of length 2.',length(K));
    end
else %K is matrix
    if isequal(size(K),[3,3])
        cam_model = 'perspective';
    else
        error('Error in depth2Normals(): K is a matrix of size %s, but it must be of size 3x3.',num2str(size(K,1)));
    end
end

if ~isempty(mask)
    
    if strcmp(cam_model, 'perspective')
        if (~exist('xx','var') || ~exist('yy','var')) || (isempty(xx) || isempty(yy))
            [xx, yy] = meshgrid( 0:size(mask,2) - 1 , 0:size(mask,1) - 1 );
            xx = xx(mask) - K(1,3);
            yy = yy(mask) - K(2,3);
        end
    elseif strcmp(cam_model, 'orthographic')
        if (~exist('xx','var') || ~exist('yy','var')) || (isempty(xx) || isempty(yy))
            xx = [];
            yy = [];
        end
    end
    
    if (~exist('Dx','var') || ~exist('Dy','var')) || (isempty(Dx) || isempty(Dy))
          if strcmp(grad_option, 'FDH')
            [Dx, Dy] = getNabla_o2p(mask, 'Forward', 'DirichletHomogeneous');
          elseif strcmp(grad_option, 'BDH')
            [Dx, Dy] = getNabla_o2p(mask, 'Backward', 'DirichletHomogeneous');
          elseif strcmp(grad_option, 'CDH')
            [Dx, Dy] = getNabla_o2p(mask, 'Central', 'DirichletHomogeneous');
          elseif strcmp(grad_option, 'FNC')
            [Dx, Dy] = getNabla_o2p(mask, 'Forward', 'NeumannConstant');
          elseif strcmp(grad_option, 'BNC')
            [Dx, Dy] = getNabla_o2p(mask, 'Backward', 'NeumannConstant');
          elseif strcmp(grad_option, 'CNC')
            [Dx, Dy] = getNabla_o2p(mask, 'Central', 'NeumannConstant');
          elseif strcmp(grad_option, 'FNH')
            [Dx, Dy] = getNabla_o2p(mask, 'Forward', 'NeumannHomogeneous');
          elseif strcmp(grad_option, 'BNH')
            [Dx, Dy] = getNabla_o2p(mask, 'Backward', 'NeumannHomogeneous');
          elseif strcmp(grad_option, 'CNH')
            [Dx, Dy] = getNabla_o2p(mask, 'Central', 'NeumannHomogeneous');
          end

    end
    
    if size(z,2)>1 %if z is not yet a vector, make it a vector
        z = z(mask);
    end
    
end

if size(z,2)>1 %if z is not yet a vector, make it a vector
    error('Error in depth2Normals\nz is assumed to be a vector if no mask is given');
end


if issparse(Dx) && issparse(Dy)
    zx = Dx*z;
    zy = Dy*z;
    % elseif size(Dx,2) == 1 && size(Dy,2) == 1%assumes that Dx and Dy are already the derivatives
    %   zx = Dx;
    %   zy = Dy;
else
    error('Error in depth2Normals\nDx and Dy should either be sparse gradient matrices or represent the derivatives of z');
end

[normals, dz, N_unnormalized] = getNormalMap(z, zx, zy, K, xx, yy, cam_model);

if nargout >= 9
    [J_n_un, J_dz] = calcJacobian(normals, z, zx, zy, K, xx, yy, Dx, Dy, cam_model);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% getNormalMap%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [N_normalized, dz, N_unnormalized] = getNormalMap(z, zx, zy, K, xx, yy, cam_model)
%variable explanation:
% z the depth as a vector
% zx = Dx*z\in Nx1, i.e. a vector;
% zy = Dy*z\in Nx1, i.e. a vector;
% K the intrinsic matrix
% xx and yy is the meshgrid as vector, where principal point is already
% taken into account, i.e. xx = xx - K(1,3) & yy = y - K(2,3)

%%
%get number of pixel in vector
pixels = size(zx,1);

% get unnormalized normals
N_unnormalized = zeros(pixels,3);
if strcmp(cam_model, 'perspective')
    N_unnormalized(:,1) = K(1,1) * zx;
    N_unnormalized(:,2) = K(2,2) * zy;
    N_unnormalized(:,3) = ( -z -xx .* zx - yy .* zy );
  %  N_unnormalized = z./(K(1,1)*K(2,2)).*N_unnormalized;
elseif strcmp(cam_model, 'orthographic')
    N_unnormalized(:,1) = zx;
    N_unnormalized(:,2) = zy;
    N_unnormalized(:,3) = -ones([pixels,1]);
else
    error('Error in depth2Normals: unknown camera model "%s"', cam_model);
end

% get normalizing constant
dz = max(eps,sqrt( sum( N_unnormalized .^ 2, 2)  ));

% normalize normals
N_normalized = bsxfun(@times, N_unnormalized, 1./dz);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calcJacobian%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J_n_un, J_dz] = calcJacobian(normals, z, zx, zy, K, xx, yy, Dx, Dy, cam_model)

% if strcmp(cam_model, 'perspective')
%   grad_dz = [ K(1,1)*normals(:,1) - xx .* normals(:,3),...%derivative of dz wrt. theta_1 (theta_1=zx)
%     K(2,2)*normals(:,2) - yy .* normals(:,3),...%derivative of dz wrt. theta_2 (theta_2=zy)
%     -normals(:,3)              ];%derivative of dz wrt. theta_3 (theta_3=z)
%
%   %some helper variables
%   f1 = K(1,1)*ones(size(xx));
%   f2 = K(2,2)*ones(size(xx));
%   c0 = zeros(size(xx));
%   c1 = ones(size(xx));
% elseif strcmp(cam_model, 'orthographic')
%   grad_dz = normals(:,1:2);%derivative of dz wrt. theta
%
%   %some helper variables
%   f1 = ones(size(normals(:,1)));
%   f2 = f1;
%   c0 = zeros(size(normals(:,1)));
%   xx = c0;
%   yy = c0;
%
% else
%   error('Error in depth2Normals: unknown camera model "%s"', cam_model);
% end

J_n_un = cell(1,size(normals,2));
npix = length(z);


 J_n_un{1} = K(1,1)*Dx;
 J_n_un{2} = K(2,2)*Dy;
 J_n_un{3} = -speye(npix) -spdiags(xx, 0, npix, npix)*Dx - spdiags(yy,0,npix,npix)*Dy;

%  J_n_un{1} = (spdiags(zx,0,npix,npix) + spdiags(z,0,npix,npix)*Dx)./K(2,2);
%  J_n_un{2} = (spdiags(zy,0,npix,npix) + spdiags(z,0,npix,npix)*Dy)./K(1,1);
%  J_n_un{3} = 1./(K(1,1)*K(2,2)).*(spdiags(-2*z - xx.*zx - yy.*zy,0,npix,npix) - spdiags(xx.*z,0,npix,npix)*Dx - spdiags(yy.*z,0,npix,npix)*Dy);

J_dz = spdiags(normals(:,1),0,npix,npix)*J_n_un{1} + spdiags(normals(:,2),0,npix,npix)*J_n_un{2} + spdiags(normals(:,3),0,npix,npix)*J_n_un{3};

end
