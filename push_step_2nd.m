function [Ez,jz,n,f,fk]=push_step_2nd(fk,kz,K,Kappaz,VZ,pos_func,phik,dpz,t,dt,z,pz,if_relativistic,externE_on,f0,nu_slow,if_apply_filter,Nfilter,nonlineardiffusioncoeff,myfilter,my_f_filter,if_E_filter,E_filter,gammaphase,tmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% push positions - first time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% push in time (semi implicit) -full step
% calculate f
fk = (fk.*pos_func);
fk = fk.*my_f_filter;
f = (ifft(ifftshift(fk,2),[],2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
if if_not_periodic
    f(:,1) = f0(:,1);
    f(:,size(f,2)) = f0(:,size(f,2));
end
%}
% density profile
%  ne1 = ne_of_x((z+t)/tmax,z/(14*pi*gammaphase)); % actually xi + t = x
%  [ne,PZ] = meshgrid(ne1,VZ(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Maxwell / Poisson equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[Ez,Elarge,jz,n]=maxwells(f,f0,invlaplacian,if_not_periodic,kz,dpz,t,dt,z,VZ,externE_on,gammaphase,tmax);
[Ez,Elarge,jz,n]=maxwells(f,f0,kz,dpz,t,dt,z,VZ,externE_on);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make everything real 
%Elarge = real(Elarge);
%f = real(f);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Apply Krook collision operator
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(abs(nu_slow))))
 %f = f - nu_slow.*(f-f0); %Unstable you donut!
% f = (f-f0).*exp(-nu_slow*dt) +f0;
%Elarge = Elarge.*exp(-nu_slow*dt) +0.0;
%Elarge = Elarge.*(1.0 - nu_slow./max(max(nu_slow)));
 end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (if_apply_filter)
    if if_apply_filter==2
        nz=max(size(z));
        fperiodic = [f(:,(nz-Nfilter+1):nz) f f(:,1:Nfilter)];
        fperiodic= conv2(fperiodic,nonlineardiffusioncoeff*myfilter,'same');
        f = fperiodic(:,(Nfilter+1):(nz+Nfilter));
    elseif if_apply_filter==1
        
        f=nonlineardiffusion(f,dt,nonlineardiffusioncoeff);
    end

end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% push in time - full momentum step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Kappa_z = -1.0*Kappaz.*(Elarge)*dt;%+VX.*(B0+Bylarge))*dt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% apply E filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (if_E_filter)
  %  fKstar = fftshift(fft(E_filter.*f,[],1),1);
  %  fKstar = fKstar.*exp(1i*Kappa_z);
  %  fKappa = fKappa +fKstar;
  f = E_filter.*f;
end

fKappa = fftshift(fft(f,[],1),1);
fKappa = fKappa.*exp(1i*Kappa_z);



fk = (ifft(ifftshift(fKappa,1),[],1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Apply Krook collision operator
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (max(max(abs(nu_slow))))
 %f = f - nu_slow.*(f-f0); %Unstable you donut!
%  ne1 = ne_of_x((z+t)/tmax,z/(14*pi*gammaphase)); % actually xi + t = x
%  [ne,PZ] = meshgrid(ne1,VZ(:,1));


 fk = (fk-f0).*exp(-nu_slow*dt) +f0;
%Elarge = Elarge.*exp(-nu_slow*dt) +0.0;
%Elarge = Elarge.*(1.0 - nu_slow./max(max(nu_slow)));
end
 
fk = fftshift(fft(fk,[],2),2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% push positions - second time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% push in time (semi implicit) -half step
% calculate f
fk = fk.*pos_func;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

