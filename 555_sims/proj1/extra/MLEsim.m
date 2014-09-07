% MLE sim
% MLE estimate of frequency, phase, amplitude (Rife & Boorstyn 74)
% K. Bell
% 11/15/99

N = 128;                  % snapshots
n = [0:1:N-1];

SNR = 10.^([[-30:5:-10] [-9:1:0] [0:5:30]]/10);
ns = length(SNR)

f = 2007;                 % actual frequency
omega = 2*pi*f;
theta = pi/6;             % actual phase

fs = 4000;                % sampling frequency
T = 1/fs;
K= 1000;                  % monte carlo trials

wse = zeros(K,ns);        % sample squared error, omega
ase = zeros(K,ns);        % sample squared error, amplitude
tse = zeros(K,ns);        % sample squared error, theta

Ct = zeros(1,ns);         % CRB, theta
Ca = zeros(1,ns);         % CRB, amplitude
Cw = zeros(1,ns);         % CRB, omega

L = 2^15;                 % FFT grid size
for m=1:ns
   alpha = sqrt(SNR(m));
   Ca(m) = 1/N;
   Cw(m) = 12/((T^2)*SNR(m)*N*(N^2-1));
   Ct(m) = 2*(2*N-1)/(SNR(m)*N*(N+1));
   for k=1:K
      x = alpha*cos(omega*n*T+theta)+randn(1,N);
      y = alpha*sin(omega*n*T+theta)+randn(1,N);
      r = x+j*y;
      
      F = fft(r,L)/N;
      P = abs(F).^2;
      [y,I] = max(P);
      omega_hat = 2*pi*fs*(I-1)/L;
      
      v = exp(j*omega_hat*n.'*T);
      f = v'*r.'/N;
      alpha_hat = abs(f);
      theta_hat = angle(f);
      
      wse(k,m) = (omega_hat-omega)^2;
      ase(k,m) = (alpha_hat-alpha)^2;
      tse(k,m) = (theta_hat-theta)^2;
   end
   [m 10*log10(SNR(m))]
end

wse_avg = sum(wse,1)/K;
ase_avg = sum(ase,1)/K;
tse_avg = sum(tse,1)/K;

figure(1)
plot(10*log10(SNR),10*log10(sqrt(wse_avg)/(2*pi)),'*')
hold on
plot(10*log10(SNR),10*log10(sqrt(Cw)/(2*pi)))
hold off
xlabel('SNR(dB)')
ylabel('RMSE (Hz)')
title('Frequency')

figure(2)
plot(10*log10(SNR),10*log10(sqrt(ase_avg)),'*')
hold on
plot(10*log10(SNR),10*log10(sqrt(Ca)))
hold off
xlabel('SNR(dB)')
ylabel('RMSE')
title('Amplitude')

figure(3)
plot(10*log10(SNR),10*log10(sqrt(tse_avg)),'*')
hold on
plot(10*log10(SNR),10*log10(sqrt(Ct)))
hold off
title('Phase')
xlabel('SNR(dB)')
ylabel('RMSE (radians)')

