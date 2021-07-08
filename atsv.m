% Completely spectral vlasov code in one dimension test code
% this version with E field solver only
% this version with filtering
% new cleaned up relativistic version 
% this version in moving frame (Correct Lorentz transform)
% This version with no normalization, correct maxwellian
function atsvrel()%
if_boosted_frame=false;
 %c = parcluster
 %matlabpool(c)
if_plotimage=true;
if_plotmovie=true;
if_relativistic = true;
if_plot_deltaf =false;
if_f_filter = false;
if_diagnose = true;

% series of filters to improve results
if_apply_filter = true; % add numerical diffusion in p/x space
numericaldiffusioncoeff = 1;%0.2;%0.05; %25; %0.05;

if_Krook=false; % slowing down of distribution - to cope with non-periodic function - function defined later
nu0 = 2.0;%0.2;

% doesn't work yet - needs conv argh!
if_E_filter = false; %stops particles accelerating from +pmax to -pmin (periodic bounds)

if_Maxwell = false; % use Maxwell solver
if_externE_on=false; % use external field

%if_boosted_frame = false; % Lorentz boost into comoving frame
gamma_boost = sqrt(1.0+1.0^2);%50.0; % lorentz factor of boost

if_dump_f_boost =false; % dump f as a funciton of time - slow and big!!!!!

system('mkdir -p movies');
system('mkdir -p movies/fld');
system('mkdir -p movies/pha');
rootdir = './movies';
% if if_boosted_frame
%     rootdir = [rootdir '2'];
% end
rootdir

% electron temperature
T=0.25;


tmax = 100;8.0*pi;50;%200.0;%100.0*(gamma_boost);

% custom time step
dtcust= 0.5;

N = 512;

nz=  N; %
npz=  N; 

ndumps=100;
kts =  (sqrt(3)/(10)^0.25/2/3.0);% sqrt(3)/2.0/sqrt(2.0)/(10)^0.25/3.0;
kbp = 1.0/(3.0/sqrt(1.0+3.0^2));

%kts = sqrt(3.0/sqrt(1.0+3.0^2))/2.0/3.0;

%kts = sqrt(3.0)/2/3;
k=1;
zmax=2.0*pi/kts;%1.0*pi/kts;
zmin=0.0;%-1.0*pi/kts;
T0=0.0001;
% max electron momentum
pzmax= 8.0;%0.6; sqrt(T0)*6;%15*sqrt(10.0);
pzmin =-8.0;% -0.6;sqrt(T0)*6;%-5.0;%-sqrt(10.0)/2.0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% wake phase velocity
v_boost = sqrt(1.0 - 1.0./gamma_boost.^2);
% boost momentum
p_boost = gamma_boost*v_boost; % 

if if_boosted_frame
    if ~if_relativistic
        message('Boost requested but not relativistic!');
    end

%  pzmax = gamma_boost*pzmax-p_boost*sqrt(1.0+pzmax^2);
 % pzmin = gamma_boost*pzmin-p_boost*sqrt(1.0+pzmin^2);
 pzmax=gamma_boost*(pzmax-v_boost*sqrt(pzmax^2+1));%3.0;%*sqrt(10);
 pzmin=gamma_boost*(pzmin-v_boost*sqrt(pzmin^2+1));%-15.0*8;
 zmax = zmax/gamma_boost;% careful with these
  zmin = zmin/gamma_boost;
end

% limits as output in picture
plim_min = pzmin;
plim_max = pzmax;
zlim_min = zmin;
zlim_max = zmax;

% distribution function
f = zeros(npz,nz);

z= linspace(zmin,zmax*(1-1/nz),nz);
pz=linspace(pzmin,pzmax,npz);

if ~if_plotimage
    filename = [rootdir '/axes'];
    save(filename, 'z','pz');
end

dz=z(2)-z(1);
dpz=pz(2)-pz(1);

[Z,PZ]=meshgrid(z,pz);

if (if_relativistic) 
    GAMMA = sqrt(1.0+PZ.^2);
else
    GAMMA = 1.0;
end


VZ=PZ./GAMMA;



%kspace
kzmax = pi/dz;
kz = linspace(-kzmax,kzmax*(1.0-2.0/nz),nz);

kappazmax = pi/dpz;

kappaz = linspace(-kappazmax,kappazmax*(1.0-2.0/npz),npz);

kappaz = kappaz;

[K,Kappaz] = meshgrid(kz,kappaz);


if if_boosted_frame
  PZlab = gamma_boost*PZ+p_boost*sqrt(1.0+PZ.^2);
  VZ = VZ + v_boost; % this shifts to the z-vpt coordinate
  tmax=tmax*gamma_boost; %careful here
  phaseboost = 0.0;%pi/2;
else
    gamma_boost=1.0;
    PZlab = PZ;
    phaseboost = 0.0;%pi;
end
%}
%initialize f
%
deltan0 = 1e-2; 1e-4;

% relativistic two stream  instability
%f0 = 1.0./sqrt(2.0*pi*T).*exp(-(abs(PZlab)-3.0).^2/2.0/T).*(1.0- deltan0.*cos(kts*gamma_boost*Z+phaseboost));
f0 = 1.0./sqrt(2.0*pi*T).*exp(-(PZ-1.5+0.1*cos(kts*Z)).^2/2/T)+1.0./sqrt(2.0*pi*T).*exp(-(PZ+1.5+0.1*cos(kts*Z)).^2/2/T);

%Landau test
%f0 = 1.0./sqrt(2.0*pi).*exp(-PZlab.^2/2.0).*(1.0- 1e-1.*cos(Z./2));


%plasma wave

%alpha=1e-1;
%f0 = 1.0./sqrt(2.0*pi*T0).*exp(-(PZlab+alpha/k*sin(Z)).^2/2.0/T0).*(1.0- alpha.*sin(Z));


% relativistic beam-plasma  instability
%f0 = 1.0./sqrt(2.0*pi*T*10).*exp(-((PZlab)-6.0*sqrt(10)).^2/2.0/T/10).*(1.0+1e-4*cos(kts*Z/gamma_boost));
%f0 = f0+1.0./sqrt(2.0*pi*T/40).*exp(-((PZlab)).^2/2.0/T*40).*(1.0+1e-4*cos(kts*Z/gamma_boost));

%f0 = f0./gamma_boost;

% This is correct shifted 1D distribution
%f0 =  gamma_phase^2*0.5/besselk(2,1.0/T,1)*(sqrt(1.0/gamma_phase^2+PZ.^2)+T/gamma_phase^2).*exp(-(gamma_phase^2/T*(sqrt(1.0/gamma_phase^2+PZ.^2)+v_phase*PZ))+1.0/T);
%f0 = neff/sqrt(2.0*pi*T)*exp(-(PZ-p_phase).^2/2.0/T);

t=0;
 text=t;

 ii=1;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % calculate E
 % calculate density
 n = sum(squeeze(sum(f0,3))*dpz,1);
n0 = n;%sum(squeeze(sum(f0,3))*dpz,1);

 f=f0;
 


 % Find dt
 dtCFL = 2*dtcust;
 if (if_Maxwell)
 dtCFL = 0.49 *dz
 end

dtwp = 0.5/sqrt(max(n))
dtcust
dt0 = min([dtCFL,dtwp,dtcust])
numericaldiffusioncoeff
nt = floor(tmax/dt0);

countdump = floor(nt/(ndumps));
counter=-1;

fKappa = zeros(npz,nz);

% this is distorted distribution function for boost
fboost = f;
if if_dump_f_boost
    fboost = zeros(npz,nz,nt);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filtering functions
%Don't exceed nz points for the filter!
% use odd number!
% this is basically numerical diffusion
Nfilter=3;

 myfilter= ones(Nfilter);
 myfilter(2,2)=-8.0;
 myfilter = myfilter*dt0/3.0 ;
 myfilter(2,2) = myfilter(2,2)+1.0;

my_f_filter = 1.0;
if if_f_filter
my_f_filter = 1.0./(1.0+exp((abs(K)/kzmax-0.9)*100.0)) ;
%my_f_filter = 1.0./(1.0+exp((sqrt((abs(K)/kzmax).^2+(abs(Kappaz)/kappazmax).^2)-0.9)*100.0)) ;
end

nu_slow=0.0; 
if if_Krook
    nu_slow = nu0*(1.0 - 1.0./(1.0+exp(((Z)-2.0/6.0*zmax)*(12*pi/zmax))+1.0*exp((-(Z)+0.9*zmin)*(12*pi/abs(zmin)))));  
   % nu_slow = nu0*(exp(-((Z-zmax).^2./2^2).^2) + exp(-((Z-zmin).^2./2^2)));
end
%plot(z,nu_slow);
%return
%E_filter = 1.0./(1.0+exp(((PZ)/pzmax-0.85)*100.0))./(1.0+exp(((PZ)/pzmin-0.85)*100.0)) -1.0;
E_filter = 1.0./(1.0+exp((PZ-pzmax+0.25)*40.0))./(1.0+exp(-(PZ-pzmin-0.25)*40.0));
%plot(pz,E_filter)
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 dnext = real(n'-1.0);%ni');

 % calculate current
 jz = sum(squeeze(sum(f.*VZ,3)),1)*dpz;

   
 nk = fftshift(fft(n-1.0));
 jzk = fftshift(fft(jz));

 %calc initial scalar potential
 phik = nk./kz.^2;
 phik = ifftshift(phik);
 phik(1) = 0.0; % zero frequency
 phik = fftshift(phik);
 phikm1 = phik;
 

 %calc electric field
 Ek = 1i*kz.*phik;

 Ek = ifftshift(Ek);
 Ez = real(ifft(Ek));


 Eext=Ez';
 [Elarge,temp] = meshgrid(Ez,pz);
 
 fk=fftshift(fft((f),[],2),2);
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %DIAGNOSTICS
 % totals on the grid
 Kene_z = sum(squeeze(sum(real(f0).*PZ.^2*0.5,3))*dpz,1);
 mass0 = sum(real(n))*dz/zmax;
 Eene0 = sum(0.5*Ez.*conj(Ez))*dz/zmax;
 Kene0 = sum(Kene_z)*dz/zmax;
 Enetot0 = Eene0+Kene0;
 
   %problem if f=0 -> add small value
   delta_ent = 1.0e-16;
 entropy0 = sum(abs(sum(squeeze(sum(f.*log(f+delta_ent),3))*dpz,1)))*dz/zmax;
 
 entropy_ext= entropy0;
 Kene_ext = Kene0;
 Eene_ext = Eene0;
 Ene_ext = Enetot0;
 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



 % need to repeat this if Maxwell on
 pos_func=exp(-1i.*K.*(VZ)*dt0*0.5);
 dumpnum=0;
 dt = dt0;

tstart = tic;
 while t<=tmax

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % get time histories of energy
 if if_diagnose  
      % totals on the grid
 Kene_z = abs(sum(squeeze(sum((f).*PZ.^2*0.5,3))*dpz,1));
 mass = sum(real(n))*dz/zmax;
 Eene = sum(0.5*Ez.*conj(Ez))*dz/zmax;
 Kene = sum(Kene_z)*dz/zmax;
 Enetot = Eene+Kene;
 
 % entropy
 %problem if f=0 -> add small value
 entropy = sum(abs(sum(squeeze(sum(f.*log(f+delta_ent),3))*dpz,1)))*dz/zmax;
 
 entropy_ext=[ entropy_ext entropy];
  Kene_ext = [Kene_ext Kene];
 Eene_ext = [Eene_ext Eene];
 Ene_ext = [Ene_ext Enetot];
 end
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% first dump output
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  counter=counter-1;
  
  

 if (counter<0)
%if (fix_positrons)
    if ~if_plot_deltaf
 pha_pz_z = real(f);
    else
         pha_pz_z = real(f)-real(f0);
    end
%end
 pha_pz_zi = flipdim(real(fi),1);
 
 if (if_plotmovie)
     npanels = 3;
   %  set(gcf, 'Colormap',...
   % [0 1 0;0.00833333376795053 0.899999976158142 0;0.0166666675359011 0.800000011920929 0;0.025000000372529 0.699999988079071 0;0.0333333350718021 0.600000023841858 0;0.0416666679084301 0.5 0;0.0500000007450581 0.400000005960464 0;0.058333333581686 0.300000011920929 0;0.0666666701436043 0.200000002980232 0;0.0750000029802322 0.100000001490116 0;0.0833333358168602 0 0;0.15384615957737 0 0;0.224358975887299 0 0;0.294871807098389 0 0;0.365384608507156 0 0;0.435897439718246 0 0;0.506410241127014 0 0;0.576923072338104 0 0;0.647435903549194 0 0;0.717948734760284 0 0;0.788461565971375 0 0;0.85897433757782 0 0;0.92948716878891 0 0;1 0 0;1 0.0416666679084301 0;1 0.0833333358168602 0;1 0.125 0;1 0.16666667163372 0;1 0.20833332836628 0;1 0.25 0;1 0.291666656732559 0;1 0.333333343267441 0;1 0.375 0;1 0.416666656732559 0;1 0.458333343267441 0;1 0.5 0;1 0.541666686534882 0;1 0.583333313465118 0;1 0.625 0;1 0.666666686534882 0;1 0.708333313465118 0;1 0.75 0;1 0.791666686534882 0;1 0.833333313465118 0;1 0.875 0;1 0.916666686534882 0;1 0.958333313465118 0;1 1 0;1 1 0.10416666418314;1 1 0.20833332836628;1 1 0.3125;1 1 0.416666656732559;1 1 0.520833313465118;1 1 0.625;1 1 0.729166686534882;1 1 0.833333313465118;1 1 0.9375;0.925490200519562 0.866106450557709 0.931302547454834;0.850980401039124 0.732212901115417 0.925105035305023;0.776470601558685 0.598319351673126 0.918907582759857;0.701960802078247 0.464425772428513 0.912710070610046;0.627451002597809 0.330532222986221 0.90651261806488;0.552941203117371 0.196638658642769 0.90031510591507;0.47843137383461 0.062745101749897 0.894117653369904]);

 subplot(npanels,1,1);
 
 imagesc(z,pz,pha_pz_z); %hold on; 
 title(['t = ' num2str(t)]);

   xlabel('z');
 ylabel('p_z');
 axis xy tight
 ylim([plim_min plim_max]); %HACK!

  subplot(npanels,1,2);
  title(['t = ' num2str(t)]);
plot(z,real(n-n0),z,real(jz));
  
axis tight
   xlabel('z');
   
legend('show')
legend('{\Delta}n','jz');
  subplot(npanels,1,3);
 plot(z,real(Ez)); 
 legend('show')
legend('E_z');
 axis tight
   xlabel('z');

     
     
 pause(0.01);
%} 
 %save phasespaces
 
 elseif ~if_plotmovie
 filename = [rootdir '/pha/phasespaces_' int2str(dumpnum)];
 save(filename, 'pha_pz_z','t','-v6');%,'pha_px_z','pha_pz_px');

 filename = [rootdir '/fld/em_fields_' int2str(dumpnum)];
 save(filename, 'Ez','t');

 filename = [rootdir '/fld/j_rho_' int2str(dumpnum)];
 save(filename, 'jz','n','t');

 
 end

 dumpnum=dumpnum+1;
counter = countdump-1;

elapsedtime = toc(tstart); %(dumpnum) = toc + elapsedtime(dumpnum-1);
disp(['output: ' num2str(dumpnum-1) ', simulation time = ' num2str(t) '/' num2str(tmax) ', elapsed time = ' num2str(elapsedtime) ' s']);
if if_diagnose
%     % totals on the grid
%  Kene_z = sum(squeeze(sum(real(f0).*PZ.^2*0.5,3))*dpz,1);
%  mass = sum(real(n))*dz/zmax;
%  Eene = sum(0.5*Ez.^2)*dz/zmax;
%  Kene = sum(Kene_z)*dz/zmax;
%  Enetot = Eene+Kene;
    disp(['dm = ' num2str(mass-mass0) ', Eene = ' num2str(Eene) ', Kene = ' num2str(Kene) ', dEne_tot = ' num2str(Enetot-Enetot0), ', dEntropy = ' num2str(entropy-entropy0)]);
end
tic;

 end
 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  
  % MAIN STEP 

%[Ez,jz,n,f,fk]=push_step_2nd(fk,kz,K,Kappaz,VZ,pos_func,phik,dpz,t,dt,z,pz,if_relativistic,if_externE_on,f0,nu_slow,if_apply_filter,Nfilter,numericaldiffusioncoeff,myfilter,my_f_filter,if_E_filter,E_filter);
[Ez,jz,n,f,fk]=push_step_2nd(fk,kz,K,Kappaz,VZ,pos_func,phik,dpz,t,dt,z,pz,if_relativistic,if_externE_on,f0,nu_slow,if_apply_filter,Nfilter,numericaldiffusioncoeff,myfilter,my_f_filter,if_E_filter,E_filter,gamma_boost,tmax);

     
 dnext = [ dnext real(n'-n0')];%ni')];
  Eext = [Eext real(Ez')];
 text(ii+1)=t;
 if if_dump_f_boost
      fboost(:,:,ii+1) = f;
 end
 
 
 ii=ii+1;
 


 t=t+dt;

 end

if if_plotimage 
    subplot(2,1,1);imagesc(text,z,dnext); colorbar; title('\deltan'); xlabel('\omega t ');ylabel('z\omega/c')
  subplot(2,1,2);imagesc(text,z,Eext); colorbar; title('E');xlabel('\omega t');ylabel('z\omega/c')
else
    filename = [rootdir '/params'];
    save(filename,'text','z','pz','Z','PZ','gamma_boost','T')
 filename = [rootdir '/fld/n_E_int'];
 save(filename, 'dnext','Eext');
filename = [rootdir '/fld/Ene_int'];
 save(filename, 'text','z','Ene_ext','Eene_ext','Kene_ext','entropy_ext');
   if if_dump_f_boost
      filename = [rootdir '/pha/fboost'];
      t=text;
 save(filename, 'fboost','t','-v6');
  end

end
 
 