% PFM 2 source CRB
% K. Bell
% 11/28/99

a = [0.01:0.01:10*pi];
e = 3*(sin(a)./a + 2*(cos(a)-sin(a)./a)./(a.^2));
a = [0 a];
e = [1 e];
figure(1)
plot(a/pi,e)
xlabel('\alpha/\pi')
ylabel('\eta(\alpha)')
grid on

E2_No = 10.^([-10:1:30]/10);
beta = 20*sqrt(3);
T = 0.5;
sa = 2*sqrt(3);
DA = [0.001 0.01 0.05 0.1 1 2]'*sa;
a = 0.5*beta*T*DA;
eta = 3*(sin(a)./a + 2*(cos(a)-sin(a)./a)./(a.^2));
CRB = ((1-eta.^2).^(-1))*((E2_No).^(-1))*(12/((beta*T)^2));
plot(10*log10(E2_No),10*log10(sqrt(CRB)))
xlabel('2E/No(dB)')
ylabel('sqrt(CRB)')