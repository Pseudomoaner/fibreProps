function [N,Lp] = im_scalablehess2(im, scale_range)
%IM_SCALABLEHESS2 generates the scale-space representation of an input
%image, finding the eigenvalues of the Hessian and the gradient pointing
%down the local axis of maximum curvature at each point.
%
%   INPUTS:
%       -im: Input image, matrix of double values
%       -scale_range: Scalar, or vector, defining range of ridge scales
%       calculation hould be performed over
%
%   OUTPUTS:
%       -N: N-score for each spatial location at each spatial scale.
%       -Lp: Component of the image gradient pointing down the local axis
%       of maximum curvature
%
%   Authors: Joeseph Harvey (c) 2014 and Oliver J. Meacock (c) 2019

gamma = 0.03; %Originally 0.05 (22/10/2020)

dim = 0;
N = zeros([size(im),size(scale_range,2)]);
Lp = zeros([size(im),size(scale_range,2)]);
for scale = scale_range
    dim = dim + 1;
    [N(:,:,dim), cosBet, sinBet] = im_hessangle2(im, scale, gamma);
    Lp(:,:,dim) = im_curvegrad2(im, scale, cosBet, sinBet);
end