function [jac] = arakawa(p,o,dx,dy)
% Function wich implements Arakawa discretisatio of Jacobian J(p,o). Fields p and o must have the same dimension.
% Grid is uniform
[M, N] = size(p);
jac=zeros(M,N);

for i=2:M-1
    for j=2:N-1
        jpp=(p(i+1,j)-p(i-1,j))*(o(i,j+1)-o(i,j-1))-(o(i+1,j)-o(i-1,j))*(p(i,j+1)-p(i,j-1));
        jpx=(p(i+1,j)*(o(i+1,j+1)-o(i+1,j-1))- p(i-1,j)*(o(i-1,j+1)-o(i-1,j-1))-p(i,j+1)*(o(i+1,j+1)-o(i-1,j+1))+p(i,j-1)*(o(i+1,j-1)-o(i-1,j-1)));
        jxp=-(o(i+1,j)*(p(i+1,j+1)-p(i+1,j-1))- o(i-1,j)*(p(i-1,j+1)-p(i-1,j-1))-o(i,j+1)*(p(i+1,j+1)-p(i-1,j+1))+o(i,j-1)*(p(i+1,j-1)-p(i-1,j-1)));
        jac(i,j)=jpp+jxp+jpx;
    end
end
jac=jac/(12*dx*dy);
end