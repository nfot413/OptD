% clear,clc
% close all;
% kr=0.6;
% kl=0.5;
% k1=1;
% k3=0.31451;
% phi1=0.5*2*pi;
% phi3=0.5*2*pi;
% phid=0;
% t=1;
% fFSR=320e9;
kr=0.5;
kl=0.5;
k1=0.8673;
k2=0.09;
k3=0.393;
phi1=0.5*2*pi;
phi2=0.5*2*pi;
phi3=0.5*2*pi;
phid=0;
t=0.979888;
fFSR=100e9;

H1=[sqrt(1-kr) -j*sqrt(kr);
   -j*sqrt(kr) sqrt(1-kr) ];
H3=[sqrt(1-kl) -j*sqrt(kl);
   -j*sqrt(kl) sqrt(1-kl) ];

Au=@(w) (sqrt(1-k1)-t^2.*exp(-j.*w).^2.*exp(-j*phi1))./(1-sqrt(1-k1).*t^2.*exp(-j.*w).^2.*exp(-j*phi1)).*(sqrt(1-k2)-t^2.*exp(-j.*w).^2.*exp(-j*phi2))./(1-sqrt(1-k2).*t^2.*exp(-j.*w).^2.*exp(-j*phi2));
Al=@(w) (sqrt(1-k3)-t^2.*exp(-j.*w).^2.*exp(-j*phi3))./(1-sqrt(1-k3).*t^2.*exp(-j.*w).^2.*exp(-j*phi3)).*t.*exp(-j.*w).^1.*exp(-j.*phid);
s=10;
w1=-20*pi;
w2=20*pi;
w=w1:0.006285:w2;
len=length(w);
w=w1-0.006285;

for i=1:len
  w=w+0.006285; 
  H2=[Au(w)  0;
    0   Al(w)];
H=H1*H2*H3;
H11(i)=H(1,1);
H12(i)=H(1,2);
H21(i)=H(2,1);
H22(i)=H(2,2);
end
frequency_f=linspace(193.1e12-s*fFSR-12.5e9,193.1e12+s*fFSR-12.5e9,len);


% figure
% plot((phase(H21)./2./pi))

hold on
amplitude_y=20*log10(abs(H11));
%  plot(frequency_f,amplitude_y,'r')
%%%%%%%%%%%%%   no phase response  %%%%%%%%%%%%%%%%%
t=length(amplitude_y);    %  know the length of the number for measured filter
% phase_y=phase(H11)./2./pi;
%%%%%%%%%%%%%%% test over %%%%%%%%%%%%%%%%
%   phase_y=phase_y.*360;


frequency_f=linspace(193.1e12-s*fFSR,193.1e12+s*fFSR,len);
% phase_y=zeros(1,19995);
plot(frequency_f.', amplitude_y)

% data_filter = [frequency_f.' amplitude_y.' phase_y.'];
% 
%  save(['VSB_RAMZI','.txt'],'data_filter','-ascii','-double');
 
%  a= fliplr(abs(fftshift(fft(H21,1000))));
%  plot(a./17.4) 