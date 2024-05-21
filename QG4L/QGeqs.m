function [q, psi, psiold] = QGeqs(n,psi,psiold,q,M,N,bcflag,bcval,...
    dissrate,BH,AH,dt,dx,dy,extra,optsur)

% From known vorticity, get streamfunction with intial guess extrapolated
    % from two previous values of psi
    if n>1
        tol=0.000001; 
        work=inverse(q,psi+extra*(psi-psiold),dx,dy,(M+N),tol,optsur,bcflag,bcval);       
        psiold=psi;
        psi=work;
    end

    % From vorticity and streamfunction, calculate jacobian
    jacobphys=arakawa(psi,q,dx,dy);


    % Laplacian of q if friction is present
    if BH~=0 || AH~=0
        lap=laplaciancr(q,dx,dy);
        if bcflag==1
            lap=periodic(lap);
        elseif bcflag==3
            lap=periodicx(lap);
        else
            lap=dirichletbc(lap,bcval);
        end

        LAP = laplaciancr(lap,dx,dy);
    else
        LAP=0;
        lap=0;
    end
    
    % Update vorticity by time evolution, with implicit dissipation
    qsp=q-dt*jacobphys -BH*dt*LAP + AH*dt*lap; 
    qsp=qsp/(1+dt*dissrate);

    % boundary conditions
    if bcflag==1
        qsp=periodic(qsp);
    elseif bcflag==2
        qsp=dirichletbc(qsp,bcval);
    else
        qsp=periodicx(qsp);
    end
    
    % For basic Euler scheme:
    %q.(ll) = qsp.(ll);

    % For RK scheme:
    %Corrector/update: same equation but with estimate new value of q 
  
    jaco=arakawa(psi,qsp,dx,dy);

    lapo=laplaciancr(qsp,dx,dy);
    if bcflag==1
    lapo=periodic(lapo);
    elseif bcflag==3
    lapo=periodicx(lapo);
    end

    LAPO = laplaciancr(lapo,dx,dy);

    q=qsp-dt*jaco -BH*dt*LAPO + AH*dt*lapo;

    q=q/(1+dt*dissrate);

    % boundary conditions
    if bcflag==1
        q=periodic(q);
    elseif bcflag==2
        q=dirichletbc(q,bcval);
    else
        q=periodicx(q);
    end
end
    