close all; clear all; clc;


%% GENERATE NOISE
d = 20;                           % duty cycle; default is 50 (half on)
snrdb = 0;
snr = 10^(snrdb/10.0);            % SNR (not in dB) 
E = 1/(d/100);%((d/100)/f).^2;    % Energy per one cycle (per one sample)
N0 = (E/snr)/4;                   % /4 is me scaling                     
sigma = sqrt(N0/2);               % Variance of noise samples


%% GENERATE SIGNAL 1
dt = 1/1e9;
f = 4e4;
num = 10;               % number of cycles
time = (1/f)*num-dt;
t = 0:dt:time;
sig1 = sin(2*pi*f*t);

fs = 2e5;               % samples per sec
n = 0:(10*fs/f-1);
sig1_s = sin(2*pi*f*n/fs);
coefs1 = sig1_s(length(sig1_s):-1:1);
signal1 = horzcat(zeros(1,100), sig1_s, zeros(1,300));
mfltr1 = filter(coefs1,1,signal1);
sq1 = mfltr1.^2;

noise=sigma*(randn(1,length(signal1)));
signal1_n = signal1 + noise;
mfltr1_n = filter(coefs1,1,signal1_n);
sq1_n = mfltr1_n.^2;


%% GENERATE SIGNAL 2
% ph = [1 0 1 0 1 1 0 0 1 1 1 0 1 0 0 0 1 1 0 0]; % corr ~ 500
% ph = [0 1 0 1 1 1 1 0 1 0 1 1 0 1 1 1 0 1 0 1]; % corr ~ 1000
ph = [1 0 0 0 0 0 0 1 0 1 0 0 0 0 1 0 1 1 0 1]; % corr ~ 1000
% ph = [0 0 0 0 1 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0 0 0 0 1 1 0 1 0 0 1 0 0]; % corr ~ 1000
% ph = [1 1 0 0 0 0 0 0 1 1 0 1 1 0 1 0 0 0 1 0 0 0 0 1 0 1 0 0 0 1 1 0]; % corr ~ 3000
fs = 2e5;                                       % samples per sec
for a = 1:2:(length(ph)-1)
    for n = 0:(fs/f-1)
        if (n == 0 && a==1)
            sig2_s = sin(2*pi*f*n/fs+pi*ph(a));
        else
            if (n < ceil(fs/f*0.5))
                sig2_s = horzcat(sig2_s,sin(2*pi*f*n/fs+pi*ph(a)));
            else
                sig2_s = horzcat(sig2_s,sin(2*pi*f*n/fs+pi*ph(a+1)));
            end
        end
    end
end

coefs2 = sig2_s(length(sig1_s):-1:1);
signal2 = horzcat(zeros(1,100), sig2_s, zeros(1,300));
mfltr2 = filter(coefs2,1,signal2);
sq2 = mfltr2.^2;

noise=sigma*(randn(1,length(signal2)));
signal2_n = signal2 + noise;
mfltr2_n = filter(coefs2,1,signal2_n);
sq2_n = mfltr2_n.^2;


%% GENERATE SIGNAL 3 (signal1 & coefs2)
mfltr3 = filter(coefs2,1,signal1);
sq3 = mfltr3.^2;

noise=sigma*(randn(1,length(signal1)));
signal1_n = signal1 + noise;
mfltr3_n = filter(coefs2,1,signal1_n);
sq3_n = mfltr3_n.^2;


%% PLOTS
% signal 1
figure(1);
subplot(4,2,1);
plot(signal1);
subplot(4,2,3);
stem(coefs1);
subplot(4,2,5);
plot(mfltr1);
subplot(4,2,7);
plot(sq1);

subplot(4,2,2);
plot(signal1_n);
subplot(4,2,4);
stem(coefs1);
subplot(4,2,6);
plot(mfltr1_n);
subplot(4,2,8);
plot(sq1_n);

% signal 2
figure(2);
subplot(4,2,1);
plot(signal2);
subplot(4,2,3);
stem(coefs2);
subplot(4,2,5);
plot(mfltr2);
subplot(4,2,7);
plot(sq2);

subplot(4,2,2);
plot(signal2_n);
subplot(4,2,4);
stem(coefs2);
subplot(4,2,6);
plot(mfltr2_n);
subplot(4,2,8);
plot(sq2_n);

% signal 3
figure(3);
subplot(4,2,1);
plot(signal1);
subplot(4,2,3);
stem(coefs2);
subplot(4,2,5);
plot(mfltr3);
subplot(4,2,7);
plot(sq3);

subplot(4,2,2);
plot(signal1_n);
subplot(4,2,4);
stem(coefs2);
subplot(4,2,6);
plot(mfltr3_n);
subplot(4,2,8);
plot(sq3_n);







%% ************************************************************

% %temp3 = 100; % how long to wait
% temp3 = 3;  % how long to wait
% wait = zeros(1,temp3*length(sig1));
% temp1 = horzcat(sig1,wait);
% disp('wait for (ms)'); disp((num/f-res)*temp3*1000);
% 
% % generate signal2
% res = 1/1e9;
% f = 4e4; %num = 10; % number of cycles
% %ph = [1 0 1 0 1 1 0 0 1 1 1 0 1 0 0 0 1 1 0 0]; % corr ~ 500
% %ph = [0 1 0 1 1 1 1 0 1 0 1 1 0 1 1 1 0 1 0 1]; % corr ~ 1000
% ph = [1 1 0 0 0 0 0 1 0 1 0 0 0 0 1 0 1 1 0 1]; % corr ~ 1000
% %ph = [0 0 0 0 1 1 0 0 1 0 1 0 0 1 1 1 1 1 0 0 0 0 0 1 1 0 1 0 0 1 0 0]; % corr ~ 1000
% % ph = [1 1 0 0 0 0 0 0 1 1 0 1 1 0 1 0 0 0 1 0 0 0 0 1 0 1 0 0 0 1 1 0]; % corr ~ 3000
% num = length(ph)/2; % # of cycles
% dt = 1/(2*f); % half cycle
% time = 0:res:dt-res;
% %t = 0:res:(dt*2*num)-res;
% t = 0:res:(dt*length(ph))-res;
% 
% sig2 = sin(2*pi*f*time+pi*ph(1));
% for b = 2:length(ph)
%     sig2 = horzcat(sig2, sin(2*pi*f*time+pi*ph(b)));
% end
% %n = 100; % how long to wait
% temp3 = 3;
% wait = zeros(1,temp3*length(sig2));
% temp2 = horzcat(sig2,wait);
% disp('wait for (ms)'); disp((num/f-res)*temp3*1000);
% 
% n=sigma*(randn(1,(temp3+1)*length(t)));
% 
% 
% if(1)   % input 1 for sig1/temp1, or input 0 for sig2/temp2
%     y = sig1; tx = temp1;
% else
%     y = sig2; tx = temp2;
% end 
% plot(t,y);
% 
% 
% 
% 
% 
% % add delay
% delay = 100000;
% rx = sqrt(E)*horzcat(zeros(1,delay),tx(1,1:end-delay)) + 0*n;
% 
% 
% 
% % sample1
% samp = 5; % # of samples per cycle
% dt = (1/f)/samp; % sample duration
% adc = zeros(1,length(y));
% index = 0;
% coefs = zeros(1,samp*num);
% for i=1:length(y)
%     if (dt*index <= t(1,i))
%         adc(1,i) = y(1,i);
%         %coefs(1,index+1) = y(1,i);
%         coefs(1,samp*num-(index-1)) = y(1,i); % time reversed
%         index = index+1;
%     else
%         adc(1,i) = adc(1,i-1); % zero-hold
%     end
% end
% 
% % sample2
% samp = 5; % # of samples per cycle
% dt = (1/f)/samp; % sample duration
% t2 = 0:res:(num/f)*(temp3+1)-res;
% rx_adc = zeros(1,length(rx));
% index = 0;
% for i=1:length(rx)
%     if (dt*index <= t2(1,i))
%         rx_adc(1,i) = rx(1,i);
%         index = index+1;
%     else
%         rx_adc(1,i) = rx_adc(1,i-1); % zero-hold
%     end
% end
% samp_rate = f*samp;
% disp('sampe rate (Sps)'); disp(samp_rate);
% 
% 
% % xcorr
% mfltr = filter(coefs,1,rx_adc);
% figure;
% subplot(3,1,1);
% plot(coefs);
% subplot(3,1,2);
% plot(rx_adc);
% subplot(3,1,3);
% plot(mfltr);
% % %%
% % %mfltr = filter(coefs,1,rx_adc);
% % mfltr = mfltr.^2;
% % 
% % %%
% % 
% % % plot
% % figure(1);
% % subplot(3,1,1);
% % plot(t,y); title('continuous rx`d signal');
% % subplot(3,1,2);
% % plot(t,adc); title('sampled rx`d signal');
% % subplot(3,1,3);
% % stem(coefs); title('coefs (time reversed)');
% % 
% % figure(2);
% % subplot(3,1,1);
% % plot(t2,rx); title('pre-sample rx`d signal');
% % subplot(3,1,2);
% % plot(rx_adc); title('rx sampled');
% % subplot(3,1,3);
% % plot(mfltr); title('matched filter^2');
% % 
% % 
% % disp('max corr'); disp(max(abs(mfltr)));
