% MPSK
% K. Bell
% 11/22/99

close all; clear; clc;



% SIMULATION PARAMETERS
M = 8;
K = log2(M);
switch M
    case 2                    % BPSK symbols
       sym_map=[1;-1];
       Es = 10.^([[-17:3:-2] [-1:1:13]]/10); % Es
    case 4                    % QPSK symbols
       sym_map=[1;1j;-1;-1j];
       Es = 10.^([[-17:3:1] [2:1:16]]/10); % Es
    case 8                    % 8PSK symbols
       sym_map=[1;(1+1j)/sqrt(2);1j;(-1+1j)/sqrt(2);-1;(-1-1j)/sqrt(2);-1j;(1-1j)/sqrt(2)];
       Es = 10.^([[-17:3:7] [8:1:22]]/10); % Es
end

Eb = Es/K;
na = length(Es);
Ns = 1000;                    % Number of symbols
No = 1;                     % noise unit variance

Es_No = Es/No;
Eb_No = Eb/No;

bits = round(rand(K,Ns));           % KxNs matrix of random 0,1 bits
symbols = graymapPSK(bits);         % fcn produces unit magnitude symbols (gray coded)



% transmit symbols over channel



% SIMULATIONS
BER = zeros(1,na); SER = zeros(1,na); Pseint = zeros(1,na);
dphi = 0.01*pi/M; phi = -pi/M+dphi/2:dphi:pi/M; nphi = length(phi); dv = 0.01;
for m=1:na
    
   r = sqrt(Es(m))*symbols + sqrt(No/2)*randn(2,Ns);       % observation eq.
   cr = r(1,:)+1j*r(2,:);                                  % converts observation eq. to complex
   sd = zeros(2,Ns);
   for n=1:Ns
      [ee ind]=min(abs(sym_map-cr(n)));                    % find minimum distance [ee=min. dist. value, ind=index of min. dist.]
      sd(:,n)= [real(sym_map(ind));imag(sym_map(ind))];    % symbol decision
   end
   bd = grayunmapPSK(sd,M);             % fcn produces bit decision
   errors = abs(bd-bits);               % differs symbol bits
   sym_err = sum(errors,1);             % accumulates # of bit errors per symbol
   sym_err = sign(sym_err);             % if 1 or more bit errors, then 1 symbol error. else 0 symbol errors
   SER(m) = sum(sym_err)/Ns;            % # of symbol errors / # of symbols
   BER(m) = sum(sum(errors))/(K*Ns);    % sums all bit errors w/in symbols, then sums up all those values for total bit error
   p1 = zeros(1,nphi);
  for q = 1:nphi
      mv = sqrt(2*Es_No(m))*cos(phi(q));
      v = max(dv/2,mv-5+dv/2):dv:mv+5;
      pv = exp(-0.5*((v-mv).^2))/sqrt(2*pi);
      p1(q) = sum(v.*pv)*dv; 
  end
   ep = exp(-Es_No(m)*sin(phi).^2).*p1/sqrt(2*pi);
   Pseint(m) = 1-sum(ep)*dphi;

end



% THEORETICAL CALCULATIONS
% erfc*(x) = 0.5*erfc(x/sqrt(2))
switch M
    case 2
%        Pse = 0.5*erfc(sqrt(2*Eb_No)/sqrt(2));
%        Pbe = Pse;
       Pbe = qfunc(sqrt(2*Eb_No));
       Pse = Pbe;
    case 4
       Pbe = qfunc(sqrt(2*Eb_No));
       Pse = 1-(1-Pbe).^2;
    case 8
       Pse = Pseint;
       Pbe = Pse/K;
end

% Pslb = 0.5*erfc(sqrt(2*Es_No)*sin(pi/M)/sqrt(2));%lower and upper bounds
% Psub = 2*Pslb;



% PLOTS
if(1)
figure(1);
subplot(1,2,1)
   semilogy(10*log10(Es_No),SER,'ok', 10*log10(Es_No),Pse,'g', ...
            10*log10(Es_No),Pseint,'or');%, ...
%             10*log10(Es_No),Pslb,'--b', 10*log10(Es_No),Psub,'--b');
   axis([-20 20 1e-5 1])
   legend('SER sim.','SER theo.','SER int.');%, ...
          %'SER LB','SER UB');
   
   title(['Symbol Error Rate, M=' int2str(M)])
   xlabel('E_s/N_o (dB)')
   ylabel('SER')
   hold off
   
subplot(1,2,2)
   semilogy(10*log10(Eb_No),BER,'ok', ...
            10*log10(Eb_No),Pbe,'-g');
   legend('BER simulated','BER theoretical');
   
   title(['Bit Error Rate, M=' int2str(M)])
   xlabel('E_b/N_o (dB)')
   ylabel('BER')
   hold off
   axis([-20 20 1e-5 1])
end   
   
   
   % EXTRA
   
   %===============================================================
   %Satellite link numerical example 
   %===============================================================
   %Pr/No=EbRb/No=EsRs/NO
   %Wc=channel bandwidth
   %Rb=bit rate
   %Pr/No=power received/noise spectral density=dBW/dBW-Hz=dB-Hz
   %Remember that Pr/No is the level of the received signal plus noise at
   %the receiver's demodulator and includes the noise figure of the
   %receiver.
   %PB=bit error rate
   %PE=probability of error
   %Rs=symbol rate
   %Es=Energy of symbol
   %Eb=Energy of bit
   %log2(4)=2
   %log2(8)=3
   %===============================================================
   %(Wc=20MHz)  (Rb=30Mbits/sec)  (Pr/No=85dB-Hz)  (PB<1e-3 required)
   %Rb>Wc therefore band-limited channel requiring MPSK scheme(QPSK or 8PSK)
   %Therefore M=4 and Rs=Rb/(log2M)=30Mbits/sec/2=15Msymbols/sec which is <Wc=20MHz
   %Es/No=(log2M)Eb/No=log2M(Pr/NoRb)=2*(10exp8.5/30e6)=21.08=13.24dB
   %Eb/No=Pr/NoRb=(10exp8.5/30e6)=10.54=10*log10(10.54)=10.23dB
   %PE(M=4)=Q(sqrt(2*Es/No)*sin(pi/M))=Q(sqrt(2*21.08)*0.3826)=Q(2.484)
   %where Q(x)=.5*erfc(sqrt(x)/sqrt(2))
   %PE(M=4)=.5*erfc(1.5786/sqrt(2))=6.5e-3
   %PB=PE(4)/log2M=6.5e-3/2=3.25e-3 !!!!PB not low enough!!!!
   %Looks like QPSK didn't do the job so lets try 8PSK
   %Therefore M=8 and Rs=Rb/(log2M)=30Mbits/sec/3=10Msymbols/sec <Wc=20MHz
   %Es/No=(log2M)Eb/No=log2M(Pr/NoRb)=3*(10exp8.5/30e6)=31.62=15dB
   %Eb/No=Pr/NoRb=(10exp8.5/30e6)=10.54=10*log10(10.54)=10.23dB
   %PE(M=8)=2*Q(sqrt(2*Es/No)*sin(pi/M))=2*Q(sqrt(2*31.62)*0.3826)=2*Q(3.042)
   %PE(M=8)=2*.5*erfc(3.042/sqrt(2))=2.4e-3
   %PB=PE(8)/log2M=2.4e-3/3=8e-4 !!!!this PB is low enough and meets the required value!!!!
   %The power at the transmitter has been increased by 15dB-13.24dB=1.76dB
   %by using 8PSK over QPSK to meet the requirements of the system.
   
  %comment out only single comments
  % function symbols = graymapPSK(bits)

%K = size(bits,1);
%N = size(bits,2);

%switch K
%case 1             % BPSK
   %% maps 0=s1, 1=s0
   %% s0 = 1, s1 = -1
   %symbols = [bits*2-1;zeros(1,N)];
   %case 2
   %% maps 00=s0, 01=s1, 11=s2, 10=s3
   %% s0 = [1;0], s1 = [0;1], s2 = [-1;0], s3 = [0;-1]
%case 3
   %% maps 000=s0, 001=s1, 011=s2, 010=s3, 110=s4, 111=s5, 101=s6, 100=s7
   %% s0 = [1;0],  s1 = 1/sqrt(2)*[1;1],  s2 = [0;1],  s3 = 1/sqrt(2)*[-1;1],
   %% s4 = [-1;0], s5 = 1/sqrt(2)*[-1;-1],s6 = [0;-1], s7 = 1/sqrt(2)*[1;-1]
   %s = sum(bits,1);
   %s_even = [1-bits(1,:)-bits(2,:);bits(2,:)-bits(1,:)];
   %s_odd = (1/sqrt(2))*[-1+2*abs(bits(3,:)-bits(1,:));-1+2*abs(bits(3,:)-bits(2,:))];
   %symbols = s_even.*([1;1]*(s==0|s==2))+s_odd.*([1;1]*(s==1|s==3));
   %end
   
   %function bits = grayunmapPSK(symbols,M)

%K = log2(M);
%N = size(symbols,2);

%switch K
%case 1             % BPSK
   %% maps 0=s1, 1=s0
   %% s0 = 1, s1 = -1
   %bits = 0.5*(symbols(1,:)+1);
   %case 2
   %% maps 00=s0, 01=s1, 11=s2, 10=s3
   %% s0 = [1;0], s1 = [0;1], s2 = [-1;0], s3 = [0;-1]
   %bits = [0.5*(1-symbols(1,:)-symbols(2,:));0.5*(1+symbols(2,:)-symbols(1,:))];
   %case 3
   %% maps 000=s0, 001=s1, 011=s2, 010=s3, 110=s4, 111=s5, 101=s6, 100=s7
   %% s0 = [1;0],  s1 = 1/sqrt(2)*[1;1],  s2 = [0;1],  s3 = 1/sqrt(2)*[-1;1],
   %% s4 = [-1;0], s5 = 1/sqrt(2)*[-1;-1],s6 = [0;-1], s7 = 1/sqrt(2)*[1;-1]
   %bits_even = [0.5*(1-symbols(1,:)-symbols(2,:));0.5*(1+symbols(2,:)-symbols(1,:));...
         %abs(symbols(2,:))];
   %bits_odd = 0.5*[1-symbols(2,:)*sqrt(2);1-symbols(1,:)*sqrt(2);sqrt(2)*abs(symbols(1,:)+symbols(2,:))];
   %s = abs(symbols(1,:))>0.7 & abs(symbols(1,:))<0.8;
   %bits = bits_even.*([1;1;1]*(~s))+bits_odd.*([1;1;1]*(s));
   %end