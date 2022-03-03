function Lp = im_curvegrad2(init_im, scale, cosBet, sinBet)

%Get spatial gradients at specified scale
g1 =  im_gradflt2(init_im, scale);

%And transform these gradients to the local coordinate space set by the
%direction of maximal curvature (beta)
Lp = g1(:,:,1).*cosBet + g1(:,:,2).*sinBet;