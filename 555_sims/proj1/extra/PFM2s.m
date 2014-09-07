% PFM 2 source
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
beta = 20*sqrt(3);       % modulation index
sigma_a = 2/sqrt(3);     % A ~ U(-sigma_a*sqrt(3),sigma_a*sqrt(3))
Amax = sigma_a*sqrt(3);
M = ceil(Amax*beta*T/pi) % # orthogonal signals

%E = 10.^([[-10:5:10] [11:1:15] [20:5:30]]/10);
E = 10.^([15:5:30]/10);
ns = length(E);
No = 2;
K= 200;              % monte carlo trials


A1se = zeros(K,ns);    % estimation error squared 
A2se = zeros(K,ns);    % estimation error squared 
Ca = zeros(1,ns);     % BCRB

deltaA = (0.01)*2*pi/(beta*T);                 % grid in A space
L = ceil(2*Amax/deltaA);                       % integer # samples
deltaA = 2*Amax/L;                             % adjust grid
AA = [-Amax+deltaA/2:deltaA:Amax-deltaA/2].';  % grid points
SA = sqrt(2/T)*sin((wc+beta*AA)*(n*Ts-T/2));   % signals

DA = sigma_a;
A1 = (rand(1,1))*Amax    % uniform >0
A2 = A1-DA
a = 0.5*beta*T*DA;
eta = 3*(sinc(a/pi) + 2*(cos(a)-sin(a/pi))/(a^2));
CRB = ((1-eta.^2).^(-1))*((E2_No).^(-1))*(12/((beta*T)^2));

for m=1:ns
   Ca(m) = (6*No)/((E(m)*T*T*beta*beta)*(1-eta^2));                    % BCRB
end
for m=1:ns
   for k=1:K
      r = sqrt(2*E(m)/T)*sin((wc+beta*A1)*(n*Ts-T/2))+...
         sqrt(2*E(m)/T)*sin((wc+beta*A2)*(n*Ts-T/2))+...
         sqrt(1/Ts)*sqrt(No/2)*randn(1,N);
      c = (r*SA.')*Ts;               % correlate
      
      % initial estimates
      [y1,I1] = max(c);
      A1_hat = AA(I1);
      c2 = c-sqrt(E(m))*sinc(0.5*beta*T*(A1_hat-AA')/pi);
      [y2,I2] = max(c2);
      A2_hat = AA(I2);
      
      % assign estimates to correct parameters
      if ((abs(A1_hat-A1)>DA/2)&(abs(A2_hat-A1)<DA/2))|((abs(A2_hat-A2)>DA/2)&(abs(A1_hat-A1)<DA/2))
         A2_hat = A1_hat;
         y2 = y1;
         A1_hat = AA(I2);
         y1 = c2(I2);
      end
      
      for q=1:4
         A1k = A1_hat; % EM
         A2k = A2_hat; % EM
         c1 = c-sqrt(E(m))*sinc(0.5*beta*T*(A2k-AA')/pi)+...
              sqrt(E(m))*sinc(0.5*beta*T*(A1k-AA')/pi);      %EM
         %c1 = c-sqrt(E(m))*sinc(0.5*beta*T*(A2_hat-AA')/pi);      %AM
         [y1,I1] = max(c1);
         A1_hat = AA(I1);
         
         c2 = c-sqrt(E(m))*sinc(0.5*beta*T*(A1k-AA')/pi)+...
              sqrt(E(m))*sinc(0.5*beta*T*(A2k-AA')/pi);      %EM
         %c2 = c-sqrt(E(m))*sinc(0.5*beta*T*(A1_hat-AA')/pi);      %AM
         [y2,I2] = max(c2);
         A2_hat = AA(I2);
      end
      
      A1se(k,m) = (A1_hat-A1)^2;
      A2se(k,m) = (A2_hat-A2)^2;
      %pause
   end
   [m 10*log10(E(m))]
end

A1se_avg = sum(A1se,1)/K;
A2se_avg = sum(A2se,1)/K;

figure(1)
plot(10*log10(2*E/No),10*log10(sqrt(A1se_avg)),'*')
hold on
plot(10*log10(2*E/No),10*log10(sqrt(A2se_avg)),'o')
plot(10*log10(2*E/No),10*log10(sqrt(Ca)))
hold off
xlabel('2E/No(dB)')
ylabel('RMSE')
title(['\sigma_a \beta T = ' num2str(sigma_a*beta*T) ', \Delta A = ' num2str(DA)])
