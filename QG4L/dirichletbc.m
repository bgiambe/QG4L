function [fpt] = dirichletbc(f,c)
% Function wich applies a constant value (c) on the borders of the field
[M, N]=size(f);

f(1,:)=c;
f(:,1)=c;
f(M,:)=c;
f(:,N)=c;
fpt=f;
end