function [lap] = laplaciancr(p,dx,dy)
% Function which implements a simple discretisation of a Laplacian of p.
% Grid must be uniform.
% Values not defined on border.
[M, N]=size(p);
lap=zeros(M,N);

lap(2:M-1, 2:N-1)=(p(3:M, 2:N-1)+p(1:M-2, 2:N-1)-2*p(2:M-1, 2:N-1))/(dx*dx)+ ...
                    (p(2:M-1, 3:N)+p(2:M-1, 1:N-2)-2*p(2:M-1, 2:N-1))/(dy*dy);
end
