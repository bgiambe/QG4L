function [fpt] = periodicx(f)
% Function which ensures periodic conditions in X direction

[M, ~]=size(f);

f(1,:)=f(M-1,:);
f(M,:)=f(2,:);
fpt=f;
end