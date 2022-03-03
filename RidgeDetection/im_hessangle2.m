function [N,cosBet,sinBet] = im_hessangle2(im, scale, gamma)
%IM_HESSANGLE2 calculates the magnitude of the eigenvalues at each location
%in an input image at the given spatial scale
%
%   INPUTS:
%       -im: Input matrix defining image, matrix of doubles
%       -scale: width (in pixels) of currently specified ridge
%       -gamma: Scale parameter (see Lindeberg 1998 Eq. 47, 48)
%   
%   OUTPUTS:
%       -A: Matrix of eigenvalue magnitudes at this spatial scale
%
%   Authors: Joeseph Harvey (c) 2014 and Oliver J. Meacock (c) 2019

%Run ridge detection filtering stages - g2 corresponds to the values of the Hessian matrix for each point in the greyscale intensity matrix.
g2 = im_hessstrflt2(im, scale);

tmp = sqrt((g2(:,:,1)- g2(:,:,3)).*(g2(:,:,1)- g2(:,:,3)) + 4 * g2(:,:,2).*g2(:,:,2));

%Components of the eigenvectors associated with the Hessian g2 (in polar coordinates)
eigvalue1 = (scale^gamma)*(g2(:,:,1) + g2(:,:,3) + tmp)/2; %Formula for eigenvalues that falls out from working through the 2x2 case (Eq. 47, 48)
eigvalue2 = (scale^gamma)*(g2(:,:,1) + g2(:,:,3) - tmp)/2;

angFrac = (g2(:,:,1) - g2(:,:,3))./tmp;
cosBet = sqrt((1+angFrac)/2); %Angle of vector at each point that points perpendicular to direction of maximal curvature
sinBet = sign(g2(:,:,2)) .* sqrt((1-angFrac)/2);

%Going to use the mtric N from Lindeberg 1998 - Eq. 49. Isolates more
%ridge-like structures specifically.
N = eigvalue1.^2 - eigvalue2.^2;
% N = eigvalue1.*(abs(eigvalue1)>=abs(eigvalue2)) + eigvalue2.*(abs(eigvalue1)<abs(eigvalue2));
N(abs(eigvalue1) < abs(eigvalue2)) = 0; %Eq. 42
N(eigvalue1 < 0) = 0; %Eq. 42