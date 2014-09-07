% Morth
%
% K. Bell
% 11/22/99


Es_No = 10.^([-10:1:20]/10);   % Es (dB)
na = length(Es_No);

Pseint = zeros(3,na);

dz = 0.01;
z = [-5:dz:5];
pz = exp(-0.5*(z.^2))/sqrt(2*pi);

% erfc*(x) = 0.5*erfc(x/sqrt(2))
M = [2 4 8];
nm = length(M);
for m=1:nm   
   K = log2(M(m));
   Eb_No(m,:) = Es_No/K;
   
   for n=1:na
      n
      p2 = (1-0.5*erfc((z+sqrt(Es_No(n)*2))/sqrt(2))).^(M(m)-1);
      Pseint(m,n) = 1-sum(p2.*pz)*dz;
   end
   
   Pb(m,:) = Pseint(m,:)*0.5*(M(m)-1)/M(m);
end

figure(1);
subplot(1,2,1)
semilogy(10*log10(Es_No),Pseint)
axis([-10 20 1e-5 1])

title(['Symbol Error Rate'])
xlabel('Es/N_o (dB)')
ylabel('SER')
hold off

subplot(1,2,2)
for m=1:nm
   semilogy(10*log10(Eb_No(m,:)),Pb(m,:))
   hold on
end
title(['Bit Error Rate'])
xlabel('Eb/N_o (dB)')
ylabel('BER')
hold off
axis([-10 20 1e-5 1])
