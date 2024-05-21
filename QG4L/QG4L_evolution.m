function [qq, pp] = QG4L_evolution(A,U, M,N,NT,epsilon,imode,...
    bcflag,bcval,dissrate,BH,AH,dt,dx,dy,L,psiscale,omegascale)

[nlayer, ~]= size(A);

% Estimate optimal overrelaxation parameter 
MM=1/(sqrt(0.5/(M*M)+0.5/(N*N)));

% boundary condition
if bcflag==1
    optsur=2/(1+2*pi/MM);
elseif bcflag==2
    optsur=2/(1+pi/MM);
else
    optsur=2/(1+1.5*pi/MM);
end

% preallocation
psi = table(zeros(M,N), zeros(M,N),zeros(M,N),zeros(M,N),...
    'VariableNames',{'Layer 1', 'Layer 2','Layer 3','Layer 4'});
q = table(zeros(M,N), zeros(M,N),zeros(M,N),zeros(M,N),...
    'VariableNames',{'Layer 1', 'Layer 2','Layer 3','Layer 4'});


% grid definition
[ic, jc] = meshgrid(1:M, 1:N);
xx = (ic-M/2)*dx/L;
yy = (jc-N/2)*dy/L;
% cylindrical coordinates
rr=sqrt(xx.^2+yy.^2); 
th=atan2(yy,xx);
% radial perturbation with 
r=rr.*(1+epsilon*cos(imode*th));

% Initial state, eddy perturbed by mode imode structure
for ll=1:nlayer
    % streamfunction   
    psi.(ll)=exp(-r.^2).*psiscale*U(ll);
    % vorticity
    q.(ll)=laplaciancr(psi.(ll),dx,dy); 
  
    if ll==1
        q.(ll)=q.(ll) + A(ll,ll)*psi.(ll).*(omegascale/psiscale) +  A(ll,ll+1)*psi.(ll+1).*(omegascale/psiscale);
    elseif ll==nlayer
        q.(ll)=q.(ll) + A(ll,ll-1)*psi.(ll-1).*(omegascale/psiscale) + A(ll,ll)*psi.(ll).*(omegascale/psiscale);
    else
        q.(ll)=q.(ll)+ A(ll,ll-1)*psi.(ll-1).*(omegascale/psiscale) + A(ll,ll)*psi.(ll) +  A(ll,ll+1)*psi.(ll+1).*(omegascale/psiscale);
    end          

    % Apply boundary conditions on fields
    if bcflag==1
        psi.(ll)=periodic(psi.(ll));
        q.(ll)=periodic(q.(ll));
    elseif bcflag==2
        psi.(ll)=dirichletbc(psi.(ll),bcval);
        q.(ll)=dirichletbc(q.(ll),bcval);
    else
      psi.(ll)=periodicx(psi.(ll));
        q.(ll)=periodicx(q.(ll));
    end
end

% preallocation for saving evolution
num_updates = 10;
ifreq=NT/num_updates;
qq = cell(1, num_updates);
pp = cell(1, num_updates);
% saving initial state
qq{1}=q;
pp{1}=psi;
k=2;

% decoupling
[B, ~]= eig(A);
Bm = inv(B);

p = psi;
o = q;
for ll=1:nlayer
    p.(ll) = Bm(ll,1)*psi.(1) + Bm(ll,2)*psi.(2) + Bm(ll,3)*psi.(3) + Bm(ll,4)*psi.(4); 
    o.(ll) = Bm(ll,1)*q.(1) + Bm(ll,2)*q.(2) + Bm(ll,3)*q.(3) + Bm(ll,4)*q.(4); 
end
psi=p;
q=o;
clear var p o;

% preallocation for psi time update
psiold=psi;
% first guess on psi for time iterations 
extra=1;


% time evolution layer-independent
for n=1:NT-1
    for ll=1:nlayer 
        [q.(ll), psi.(ll), psiold.(ll)] = QGeqs(n,psi.(ll),psiold.(ll),q.(ll),M,N,bcflag,bcval,dissrate,BH,AH,dt,dx,dy,extra,optsur);
    end

  % re-coupling
    psiphys = psi;
    qphys = q;
    for ll=1:nlayer
        psiphys.(ll) = B(ll,1)*psi.(1) + B(ll,2)*psi.(2) + B(ll,3)*psi.(3) + B(ll,4)*psi.(4); 
        qphys.(ll) = B(ll,1)*q.(1) + B(ll,2)*q.(2) + B(ll,3)*q.(3) + B(ll,4)*q.(4);
    end
    
    % saving for evolution

    if mod(n,ifreq)==0
        qq{k}=qphys;
        pp{k} = psiphys;
        k=k+1;
    end
end

end