% PFM
% K. Bell
% 10/28/99


fc = 100;            % Carrier Frequency (Hz)
wc= 2*pi*fc;

fs = 1000;           % Sampling Frequency (Hz)
Ts = 1/fs;           % Sampling period (sec)

T = 0.5;             % Observation interval (sec.) >> 1/wc
N = round(T/Ts);     % Samples 
T = N*Ts;            % corrected T
n = [0.5:1:N-0.5];   % n+0.5

                         % |beta*Amax| << wc
beta = 4*sqrt(3);       % modulation index
sigma_a = 2/sqrt(3);     % A ~ U(-sigma_a*sqrt(3),sigma_a*sqrt(3))
Amax = sigma_a*sqrt(3);
M = ceil(Amax*beta*T/pi) % # orthogonal signals

E = 10.^([[-10:5:5] [6:1:16] [20:5:30]]/10);
ns = length(E);
No = 2;
K= 5000;              % monte carlo trials


Ase = zeros(K,ns);    % estimation error squared 
Ca = zeros(1,ns);     % BCRB
mse = zeros(1,ns);    % mse approx.
po = zeros(1,ns);     % po

deltaA = (0.01)*2*pi/(beta*T);                 % grid in A space
L = ceil(2*Amax/deltaA);                       % integer # samples
deltaA = 2*Amax/L;                             % adjust grid
AA = [-Amax+deltaA/2:deltaA:Amax-deltaA/2].';  % grid points
SA = sqrt(2/T)*sin((wc+beta*AA)*(n*Ts-T/2));   % signals


dz = 0.01;                                     % numerical integration for p0
z = [-5:dz:5];
pz = exp(-0.5*(z.^2))/sqrt(2*pi);

for m=1:ns
   Ca(m) = (6*No)/(E(m)*T*T*beta*beta);                    % BCRB
   p2 = (1-0.5*erfc((z+sqrt(E(m)*2/No))/sqrt(2))).^(M-1);  % po calc  
   po(m) = 1-sum(p2.*pz)*dz;
   mse(m) = 2*po(m)*sigma_a^2+ (1-po(m))*Ca(m);            % approx mse
end
for m=1:ns
   for k=1:K
      A = (rand(1,1)-0.5)*2*Amax;    % uniform
      r = sqrt(2*E(m)/T)*sin((wc+beta*A)*(n*Ts-T/2))+sqrt(1/Ts)*sqrt(No/2)*randn(1,N);
      c = (r*SA.')*Ts;               % correlate
      [y,I] = max(c);
      A_hat = AA(I);
      Ase(k,m) = (A_hat-A)^2;
   end
   [m 10*log10(E(m))]
end

Ase_avg = sum(Ase,1)/K;

figure(1)
plot(10*log10(2*E/No),10*log10(sqrt(Ase_avg)),'*')
hold on
plot(10*log10(2*E/No),10*log10(sqrt(Ca)))
plot(10*log10(2*E/No),10*log10(sqrt(mse)),'g')
hold off
xlabel('2E/No(dB)')
ylabel('RMSE')
title(['\sigma_a \beta T = ' num2str(sigma_a*beta*T)])
