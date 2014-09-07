%JC 6/20/06
%I was interested in the theory for 8PSK and came across this m-file written in 1999. It
%requires 2 function calls (graymapPSK and grayunmapPSK) which I have
%included and must be uncommented(single comments only) and loaded into your workspace.
%Remove the function calls from this m-file after you load into your
%workspace.
%It seems to work well and give valid results for BPSK and QPSK but I have a question
%about the setup for 8PSK and the solution for Pseint where Pse=Pseint.
%Also, the Gray coding seems different than what I'm used to.
%Some of you more theoretical minded may have some insight into the theory
%behind how Pseint was determined and comment. I don't have the
%communications toolbox-blockset and was wondering how it compared with Matlab's
%example simulation of 8PSK with Gray coding (in Simulink?)
%Each run takes about a minute and the loop iterations are shown in the
%command window.
%A  numerical example of a satellite link is shown using uncoded QPSK or 8PSK when the Rb is
%greater than the channel bandwidth Wc (Band-limited channel). It's a
%little wordy but I have tryed to provide enough information so you don't
%have to read between the lines.
% MPSK
% K. Bell
% 11/22/99
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 8;
K = log2(M);
switch M
case 2                    % BPSK symbols
   sym_map=[1;-1];
   Es = 10.^([[-17:3:-2] [-1:1:13]]/10); % Es
case 4                    % QPSK symbols
   sym_map=[1;j;-1;-j];
   Es = 10.^([[-17:3:1] [2:1:16]]/10); % Es
case 8                    % 8PSK symbols
   sym_map=[1;(1+j)/sqrt(2);j;(-1+j)/sqrt(2);-1;(-1-j)/sqrt(2);-j;(1-j)/sqrt(2)];
   Es = 10.^([[-17:3:7] [8:1:22]]/10); % Es
end

Eb = Es/K;
na = length(Es);
Ns = 100000;              % Number of symbols
No = 2;                   % noise unit variance

Es_No = Es/No;
Eb_No = Eb/No;

bits = round(rand(K,Ns));           % KxNs matrix of random 0,1 bits
symbols = graymapPSK(bits);         % gray code map to symbols
BER = zeros(1,na);
SER = zeros(1,na);
Pseint = zeros(1,na);

dphi = 0.01*pi/M;
phi = [-pi/M+dphi/2:dphi:pi/M];
nphi = length(phi);
dv = 0.01;

for m=1:na
   m
   r = sqrt(Es(m))*symbols + sqrt(No/2)*randn(2,Ns);         % observations
   cr = r(1,:)+j*r(2,:);
   sd = zeros(2,Ns);
   for n=1:Ns
      [ee ind]=min(abs(sym_map-cr(n)));                    % map to closest symbol
      sd(:,n)= [real(sym_map(ind));imag(sym_map(ind))];    % symbol decision
   end
   bd=symbols;
   bd = grayunmapPSK(sd,M);
   errors = abs(bd-bits);                % KxNs matrix of bit errors
   symb_err = sign(sum(errors,1));       % one symbol error per column
   SER(m) = sum(symb_err)/Ns;
   BER(m) = sum(sum(errors))/(K*Ns);
   p1 = zeros(1,nphi);
  for q = 1:nphi
      mv = sqrt(2*Es_No(m))*cos(phi(q));
      v = [max(dv/2,mv-5+dv/2):dv:mv+5];
      pv = exp(-0.5*((v-mv).^2))/sqrt(2*pi);
      p1(q) = sum(v.*pv)*dv; 
  end
   ep = exp(-Es_No(m)*sin(phi).^2).*p1/sqrt(2*pi);
   Pseint(m) = 1-sum(ep)*dphi;

end

% erfc*(x) = 0.5*erfc(x/sqrt(2))

switch M
case 2
   
   Pse = 0.5*erfc(sqrt(2*Eb_No)/sqrt(2));
   Pbe = Pse;

case 4
   Pbe = 0.5*erfc(sqrt(2*Eb_No)/sqrt(2));
   Pse = 1-(1-Pbe).^2;

case 8
   Pse = Pseint;
   Pbe = Pse/K;
end

Pslb = 0.5*erfc(sqrt(2*Es_No)*sin(pi/M)/sqrt(2));%lower and upper bounds
Psub = 2*Pslb;


figure(1);
subplot(1,2,1)
   semilogy(10*log10(Es_No),SER,'oc')
   hold on
   semilogy(10*log10(Es_No),Pse,'g')
   semilogy(10*log10(Es_No),Pseint,'b')
   semilogy(10*log10(Es_No),Pslb,'--r')
   semilogy(10*log10(Es_No),Psub,'--r')
   axis([-20 20 1e-5 1])
   
   title(['Symbol Error Rate, M=' int2str(M)])
   xlabel('Es/N_o (dB)')
   ylabel('SER')
   hold off
   
subplot(1,2,2)
   semilogy(10*log10(Eb_No),BER,'oc')
   hold on
   semilogy(10*log10(Eb_No),Pbe,'-g')
   
   title(['Bit Error Rate, M=' int2str(M)])
   xlabel('Eb/N_o (dB)')
   ylabel('BER')
   hold off
   axis([-20 20 1e-5 1])
   
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