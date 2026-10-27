function [x, w] = gauss_nodes(n)
% GAUSS_NODES  Gauss-Legendre nodes and weights on [-1,1] (lookup table
%   from the lecture, slide 15). n+1 nodes are returned, n = 0,...,5.
%   Exact for polynomials of degree <= 2n+1.
    if(n==0)
        x = 0;   w = 2;
    elseif(n==1)
        x(1) = -1/sqrt(3);      x(2) = -x(1);
        w(1) =  1;              w(2) =  w(1);
    elseif(n==2)
        x(1) = -sqrt(3/5);      x(2) = 0;           x(3) = -x(1);
        w(1) =  5/9;            w(2) = 8/9;         w(3) =  w(1);
    elseif(n==3)
        x(1) = -0.861136311594053;     x(4) = -x(1);
        x(2) = -0.339981043584856;     x(3) = -x(2);
        w(1) =  0.347854845137454;     w(4) =  w(1);
        w(2) =  0.652145154862546;     w(3) =  w(2);
    elseif(n==4)
        x(1) = -0.906179845938664;     x(5) = -x(1);
        x(2) = -0.538469310105683;     x(4) = -x(2);
        x(3) =  0;
        w(1) =  0.236926885056189;     w(5) =  w(1);
        w(2) =  0.478628670499366;     w(4) =  w(2);
        w(3) =  0.568888888888889;
    elseif(n==5)
        x(1) = -0.932469514203152;     x(6) = -x(1);
        x(2) = -0.661209386466265;     x(5) = -x(2);
        x(3) = -0.238619186083197;     x(4) = -x(3);
        w(1) =  0.171324492379170;     w(6) =  w(1);
        w(2) =  0.360761573048139;     w(5) =  w(2);
        w(3) =  0.467913934572691;     w(4) =  w(3);
    else
        error('gauss_nodes: only n = 0,...,5 are tabulated.');
    end
end
