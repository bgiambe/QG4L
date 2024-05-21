function [psif] = inverse(om,psistart,dx,dy,nmax,tol,optsur,bcflag,bcval)
% Inversion of equation Laplacian psi  = omega
% starting field for iterations is psistart
% Grid spacing is uniform dx and dy
% Maximum number of iterations is nmax
% tol is the relative error below which iterations are stopped
% optsur: overrelaxation parameter (between 0 and 2)
% bcflag: boundary flag: 1 for perdiocicity on both direction
%                        2 for dirichlet condition with value bcval
%                        3 for periodicity in x direction and fixed value
%                        in y boundaries

[M, N] = size(om);

% Initial value
psi=psistart;

% Adapt stopping criteria
n=0;
relerr=1;
err=zeros(M,N);

while n < nmax && relerr > tol
    errtot=0;
    psitot=0;
    avrpsi=0;
   
    inv_dx2=1/(dx*dx);
    inv_dy2=1/(dy*dy);
    denom=2*(inv_dx2+inv_dy2);
    
    for i=2:M-1
        for j=2:N-1
            dx2=(psi(i+1,j)+psi(i-1,j)-2*psi(i,j))*inv_dx2;
            dy2=(psi(i,j+1)+psi(i,j-1)-2*psi(i,j))*inv_dy2;
            err(i,j)=(om(i,j)-(dx2+dy2))/denom;
            % Gauss Seidel approach update
            psi(i,j)=psi(i,j)-optsur*err(i,j);
            errtot=errtot+err(i,j)*err(i,j);
            psitot=psitot+psi(i,j)*psi(i,j);
            avrpsi=avrpsi+psi(i,j);
        end
    end
   
    errtot=errtot/(N-2)/(M-2);
    psitot=psitot/(N-2)/(M-2)-(avrpsi/(N-2)/(M-2))^2;
    % For stopping criterium
    relerr=sqrt(errtot/psitot);
   
    if bcflag==1
        psi=periodic(psi);
    elseif bcflag==2
        psi=dirichletbc(psi,bcval);
    else 
    psi=periodicx(psi);

    psi(:,1)=psistart(:,1);
    psi(:,N)=psistart(:,N);
    end
    
    n=n+1;
    
    
end

psif=psi;

end