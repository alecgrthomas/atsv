function [Ez,Elarge,jz,n]=maxwells(f,f0,kz,dpz,t,dt,z,VZ,externE_on,gammaphase)

% use spectral method or FD
spectral = true;

% calculate E
 % calculate density
 n = sum(squeeze(sum(f,3))*dpz,1);
  n0 = sum(squeeze(sum(f0,3))*dpz,1);
 
   % and current
   jz = sum(squeeze(sum(f.*VZ,3)),1)*dpz;
     
 if spectral  
   
 nk = fftshift(fft(n));
 jzk = fftshift(fft(jz));

 %calc initial scalar potential
 phik = -nk./kz.^2;
 phik = ifftshift(phik);
 phik(1) = 0.0; % zero frequency offset
 phik = fftshift(phik);

 %calc electric field
 Ek = 1i*kz.*phik;
 Ek = ifftshift(Ek);
 Ez = (ifft(Ek));
 
 else % FD
    % awful quickly implemented solution - fix later
    Nz = length(kz);
    dz = pi/max(abs(kz));
    A = zeros(Nz,Nz);
    A(1,1) = -2; %open
    A(1,2) = 1;
    for ii = 2:Nz-1;
    A(ii,ii) = -2;
    A(ii,ii+1) = 1;
    A(ii,ii-1) = 1;
    end
    A(Nz,Nz-1) = 1;
    A(Nz,Nz) = -2;
    
    A=A/dz^2;
    A = inv(A);
    phi = A*(n'-mean(n0));
    fminus = circshift(phi,1.0*[1 0]);
    fplus = circshift(phi,-1.0*[1 0]);
    fminus(1,1) = -1.0*phi(1,1);
    fminus(Nz,1) = -1.0*phi(Nz,1);
    Ez = ((fplus-fminus)/(2.0*dz))';
 end 
 
% add external electric field ala Afeyan
if(externE_on)
 Ez=Ez+externE(t,z,gammaphase);
end

[Elarge,PZ] = meshgrid(Ez,VZ(:,1));



