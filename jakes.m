%%% PRABHAT KUMAR RAI --- EE18MTECH01005
%%% Implementation of Jakes Model for fc=800 MHz and v=30 Kmph

clc; close all; clear all;
N = 50; Cn = 1/sqrt(N); E = 1; T = 0.001; Ts = 5;
wm = 140; % beta*v; 50 in hertz; 2*pi*f/c *v=139.6
phi = -pi + 2*pi.*rand(1,N);
theta = -pi + 2*pi.*rand(1,N);
JakeCommat=[];
for t = 1 : T : Ts
    JakeRe = 0; JakeIm = 0; 
    for k = 1 : N
        alpha(k) = (2*pi*k-pi+theta(k))/( 4*N + 2 );
        Xc = E * Cn * (cos(wm*t*cos(alpha(k)) + phi(k)));
        Xs = E * Cn * (sin(wm*t*cos(alpha(k)) + phi(k)));
        JakeRe = JakeRe + Xc;
        JakeIm = JakeIm + Xs;
    end    
    JakeCom = JakeRe + 1i*JakeIm;
    JakeCommat = [JakeCommat JakeCom];
end

% AutoCorrelation 
real_jake = real(JakeCommat);
imag_jake = imag(JakeCommat);
envelope_jake = sqrt(real_jake.^2 + imag_jake.^2);
[a1,a2] = xcorr(JakeCommat,'Coeff');
zayz = xcorr(JakeCommat,'Coeff');

% Power Spectral Density
psdx = [-8000:8000];
Fres = 1000; N = length(zayz); psdd = fft(zayz);
xdft = psdd(1:N/2+1);
psdx = (1/(Fres*N)) * abs(xdft).^2;
% psdx(-200:-1) = 1;

psdx(2:end-1) = 2*psdx(2:end-1);
cc2 = psdx(200:-1:1);
psdx((end-199):end) = cc2;
freq = 0:Fres/length(psdd):Fres/2;

% Plotting of Autocorrelation & Power Spectral Density
len = length(a2);
bb = (len-2001)./2;
cc = bb+1 : 1 : bb+2001;
dd = -1000 : 1 : 1000;
tdd = dd*T; 
z1 = wm*tdd;  sigma0 = 1;
Jakes_bessel = sigma0.^2.*besselj(0,z1);
figure;
subplot(2,1,1);
% plot(tdd,real(a1(cc)),'-',tdd,real(T_bessel),'.');
plot(real(a1(cc)),'-','Linewidth',2); hold on;
plot(real(Jakes_bessel),'.','Linewidth',2);
title('AutoCorrelation of Jakes Model'); xlabel('x-axis not scaled')
legend('simulated','theoritical');
axis([0 2000 -0.8 1.2]); grid on;

subplot(2,1,2); 
plot(abs(psdd))
% plot(freq,psdx); grid on;
title('PSD of AutoCorrelation of Jakes Model')
xlabel('Frequency (Hz)')
axis([0 400 0 500]); grid on;



% % % 
% % % %
% % % co1 = 1 : 1000;
% % % subplot(2,2,3);
% % % semilogy(co1*T,envelope_jake(1:1000)); 
% % % % semilogy(envelope_jake); 
% % % title('Rayleigh Coefficient'); xlabel('t(second)');
% % % length_r = length(envelope_jake);
% % % pdf_env = zeros(1,501);
% % % count=0;
% % % temp = round(100.*envelope_jake);
% % % for k = 1 : length_r
% % %     if temp(k) <= 500
% % %         count = count + 1;
% % %         pdf_env(1,temp(k)+1) = pdf_env(1,temp(k)+1)+1;
% % %     end
% % % end
% % % pdf_env = pdf_env./count./0.01;
% % % 
% % % sgma2 = 0.5;
% % % x = [0:0.01:5];
% % % pdf_theory = (x./sgma2).*exp(-1.*x.^2./(2.*sgma2));
% % % 
% % % subplot(2,2,4);
% % % plot(x,pdf_env,'-',x,pdf_theory,'*');
% % % legend('Simulated','Theory'); title('PDF of jakes envelope');
