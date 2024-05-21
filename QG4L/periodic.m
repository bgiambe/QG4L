function [fpt] = periodic(f)
% Function which ensures periodic conditions

[M, N]=size(f);

f(1,:)=f(M-1,:);
f(:,1)=f(:,N-1);
f(M,:)=f(2,:);
f(:,N)=f(:,2);
fpt=f;

end