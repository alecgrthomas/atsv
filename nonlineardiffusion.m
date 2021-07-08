%apply diffusion filter then use result to feed back into new filter with
%nonlinear correction to limit gradient smoothing
function NLD=nonlineardiffusion(f,dt0,numericaldiffusioncoeff)
%numericaldiffusioncoeff = 0.05;%025;
%dt0 =0.025;
% find smoothed gradient
% 
Nfilter=3;
 myfilter= ones(Nfilter);
 myfilter(2,2)=-8.0;
 myfilter = myfilter/3.0 ;
 myfilter(2,2) = myfilter(2,2)+1.0;
 
     nz=max(size(f,2));
     fperiodic = [f(:,(nz-Nfilter+1):nz) f f(:,1:Nfilter)];
     fperiodic= conv2(fperiodic,myfilter,'same');
     u = fperiodic(:,(Nfilter+1):(nz+Nfilter));
     
% now obtain |grad u|^2

gradux = diffn(u,1,[0 1]);
graduy = diffn(u,1,[1 0]);

g = 1.0./(1.0+ (gradux.^2+graduy.^2)./(u.^2+1e-16));

gip = 0.5*(g+circshift(g,[0 1]));
gim = 0.5*(g+circshift(g,[0 -1]));
gjp = 0.5*(g+circshift(g,[1 0]));
gjm = 0.5*(g+circshift(g,[-1 0]));
%{
gip = 1;
gim =1;
gjp =1;
gjm = 1;
%}

 NLD=-1.0*(gip+gim+gjp+gjm).*f;
 NLD= NLD+gip.*circshift(f,[0 1]);
 NLD= NLD+gim.*circshift(f,[0 -1]);
 NLD= NLD+gjp.*circshift(f,[1 0]);
 NLD= NLD+gjm.*circshift(f,[-1 0]);
 
 NLD = f+NLD*0.25*numericaldiffusioncoeff*dt0;
