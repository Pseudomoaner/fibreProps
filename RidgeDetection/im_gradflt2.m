function g1 =  im_gradflt2(init_im, scale)
%IM_HESSSTRFLT2 calculates the gradient vector of the given input image, at
%the given spatial scale
%
%   INPUTS:
%       -init_im: The input image, a matrix of doubles
%       -scale: The currently specified scale over which smoothing should
%       be apllied
%
%   OUTPUTS:
%       -g1: The components of the gradient vector at each spatial location in the
%       image. First two dimensions correspond to the dimensions of the
%       input image, while three components in the third dimension are the
%       components of the Hessian (noting its symmetry).
%
%   Authors: Oliver J. Meacock (c) 2020

%Create normalized kernals - x^2, y^2 and xy (horizontal, vertical and diagonal edge detectors)
[y,x]=meshgrid(-3:6/scale:3,-3:6/scale:3);
g1a=-x.*exp(-(x.^2+y.^2));
g1a=g1a/sum(abs(g1a(:)));

g1b=-y.*exp(-(x.^2+y.^2));
g1b=g1b/sum(abs(g1b(:)));

%Filter image with edge detectors
g1a_rst=imfilter(init_im, g1a, 'symmetric', 'same');
g1(:,:,1)=g1a_rst;
g1b_rst=imfilter(init_im, g1b, 'symmetric', 'same');
g1(:,:,2)=g1b_rst;
