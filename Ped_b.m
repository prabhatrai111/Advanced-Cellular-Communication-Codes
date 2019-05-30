% PRABHAT KUMAR RAI --- EE18MTECH01005

clc; clear all; close all;

% Pedestrian-B channel
Bw = [30, 200, 20000]; % In KHz bandwidth at the reciver
P_db = [0, -0.9, -4.9, -8.0, -7.8, -23.9]; % Power in dB scale
amp = sqrt( 10.^( 0.1*P_db ) );

% Taking the sampling duration as 5ns
Td = 5; % nano second at the DSP
tds = [0, 200, 800, 1200, 2300, 3700]; % in nano sec
t = tds./Td;
h( t + 1 ) = amp;
H = fft(h);

% window function in frequency domain
z1 = 1/(Td*10^-6.*Bw(1));
z2 = 1/(Td*10^-6.*Bw(2));
z3 = 1/(Td*10^-6.*Bw(3));
z = [ z1 z2 z3 ];
p = [];
for i = 1:length(Bw)
        L = ceil( length(h)/(10*z(i)) ); % for bandwidth at channel
        W = [ ones( 1, L ), zeros(1, ( length(H) - L ) )];% Window function
        Rx = H.*W;
        rx = ifft( Rx );
        hx = abs( rx );
        h_rx = hx( 1: 100: end ); 
        subplot(3,3,i)
        stem(h_rx); xlabel( 'Number of taps' );
        grid on; ylabel( 'Amplitude' );        
        title( sprintf('Channel delay for %d KHz',Bw(i)));
   
        
        % input signal generation with 0 or 1'
        N = 10000; tap = 10;
        x_in = rand(1,N) > 0.5; 
        % bpsk generation 
        x_bpsk = 2 * x_in - 1; % conversion into bpsk bits i.e. (1 or -1)
        % upsampling the bpsk signal
        samp_fact = 5;
        x_up = upsample( x_bpsk, samp_fact );  
        
        % sinc pulse shaping
        % ray = 1/sqrt( 2*tap )*( randn( 1, tap ) + 1i*randn( 1, tap ) );
        ray = h_rx; 
        lin = -2:0.5:2;
        sinc_y1 = sinc( lin );
        sinc_y  = sinc_y1 / norm( sinc_y1 );
        new_y1 = conv( ray, sinc_y );
        new_y = new_y1/norm(new_y1);
        sinc_ray_shap = conv( x_up, new_y ); % sinc pulse shaping
        
        %-----NOISE GENERATION & ADDITION TO SIGNAl
        noise = randn( 1, length( sinc_ray_shap ) );
        comp_noise = 1/sqrt(2)*( noise + 1i*noise ); % guassian noise mean=0 variance=1 (approximately)
        Eb_No_dB = [ 0 : 20 ]; % different Eb/N values
        
        for r = 1:length( Eb_No_dB )
            
            % Noise addition
            corrp = sinc_ray_shap + 10^( -Eb_No_dB(r)/20 )*comp_noise;
            
            %----Matched filter
            match = conj( fliplr( new_y ));
            h_match = conv( match, corrp );
            
            % choosing only best samples of length l*N
            h_match = real( h_match );
            h_best = h_match( ( ( length( h_match )/2 + 1 )-( samp_fact*N/2 ) ):( ( length( h_match )/2 ) + ( samp_fact*N/2 ) ) );  
            
            % down sampling
            k_down = downsample( h_best, samp_fact );   % down sampling by 2
            
            % Demodulating the signal
            x_out = k_down > 0;
            
            % calculating the number of bit errors
            nError( r ) = biterr( x_in , x_out );
            simulation_Ber = nError/N; % simulated ber
        end
        subplot(3,3,i+3)
        semilogy( Eb_No_dB, simulation_Ber, '.-' ); hold on;
        theory_Ber = 0.5*erfc( sqrt( 10.^( Eb_No_dB/10 ) ) ); % theoretical ber 
%       theory_Ber = 0.5*( 1 - sqrt( [10.^( Eb_No_dB/10 )]/[2 + (10.^( Eb_No_dB/10 ))]));
        semilogy( Eb_No_dB, theory_Ber, 'b*-' ); hold on;
        legend( 'Simulation', 'Theory' ); xlabel( 'Eb/No, dB' );
        axis( [0 20 10^-5 0.5] ); grid on; ylabel( 'Bit Error Rate' ); 
        title( sprintf('BER curve for BPSK pulse shaping for %d KHz',Bw(i)));
       
end

        
        
        