%% PRABHAT KUMAR RAI --- EE18MTECH01005  %%

clc; clear all; close all;

%% BPSK modulation
N = 1000000;
input = rand(1, N) > 0.5; % generating 0,1 with equal probability
input_bpsk = 2*input - 1; % BPSK modulation 0 = -1, 1 = 1
samp_fact = 2;  % upsampling factor
up_input_bpsk = upsample(input_bpsk, samp_fact);

%% raised cosine pulse
T = 1; q = 5; alpha = 0.2; % Roll_off
% T = 1; q = 1; alpha = 0.2; % Roll_off
t = -q*T : 1/samp_fact : q*T;
rc_pulse = (rc_cos(alpha, t))/norm(rc_cos(alpha, t)); 
% rc_pulse = rc_cos(alpha, t);

%% Convolution of upsampled symbol and Pulse_shaping filter 
channel = [1 1];
tk = mod(length(channel), 2);
Pulse_channel = conv(channel, rc_pulse); % pulse_shape + channel
input_chan_pulse = conv(Pulse_channel, up_input_bpsk); % Convolution of upsampled symbol and Pulse_shaping filter 

%% Formation of P Toeplitz Matrix size(Nf*samp_fact, Nf+v)
p1 = []; P_toepl = []; P_toepll = []; W_MMSE_Eq = []; Z_MMSE = [];
P_even = Pulse_channel(1:2:end);
P_odd = [Pulse_channel(2:2:end) zeros(1,tk)];
p1 = [P_even; P_odd];        
v = length(Pulse_channel)/samp_fact;
Nf = 50;
for kk = 0 : (Nf - 1)
    P_toepll = [zeros(2, kk) p1 zeros(2, Nf-kk)];
    P_toepl = [P_toepl; P_toepll]; % Toeplitz matrix of Channel
end

%% BER Calculation
Eb_N0_dB = 0 : 15;

for ll = 1 : length(Eb_N0_dB)
    
    %% Noise Addition
    SNR = 10^(Eb_N0_dB(ll)/10);
    noise = sqrt(1/(2*SNR))*complex(randn(1,length(input_chan_pulse)), randn(1,length(input_chan_pulse)));
    Y_rcvd = input_chan_pulse + noise;
   
    %% MMSE equalizer
    delta = ceil((Nf+v)/2);
    delta_1 = [zeros(1, delta) 1 zeros(1, (Nf + v - delta - 1))];
    Ryy = (P_toepl'*P_toepl) + samp_fact*(1/SNR).*eye(Nf+v);
    W_MMSE_Eq = delta_1*inv(Ryy)*(P_toepl)';
   
    %% Data passing with MMSE equalizer
    Z_MMSE = conv(W_MMSE_Eq, Y_rcvd);
    Zk_even = Z_MMSE(1 : 2 : end);
    Zk_odd = Z_MMSE(2 : 2 : end);
    Zk_mat = Zk_even + Zk_odd;
    for l = 1 : N
        recvd(l) = Zk_mat(delta + l);
    end
    
    %% Output Hard Decoding
    output = real(recvd) > 0;
    BER(ll) = biterr(input, output)/N;
    Theory_BER(ll) = qfunc(sqrt(2*SNR));
end

%% figure
semilogy(Eb_N0_dB, BER, 'bs-','Linewidth',1.5); hold on;
semilogy(Eb_N0_dB, Theory_BER, 'r-*','Linewidth',2);
axis([Eb_N0_dB(1) Eb_N0_dB(end) 10^-5 0.5]); grid on;
legend('simulation mmse', 'theory'); xlabel('Eb/No, dB');
ylabel('Bit Error Rate'); title('BER curve for BPSK with MMSE Equalizer');
