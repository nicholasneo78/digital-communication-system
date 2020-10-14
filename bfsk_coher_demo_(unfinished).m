
%Define carrier frequency in Hz
Fc = 10000; 
%Given 16 times over-sampled Fs = Fc * 16
Fs = Fc * 16;
%Data rate 1kpbs
dataRate = 1000;
%Signal length
numOfBits = 1024;

%Define amp
amplitude = 1;

%low pass butter filter, 6th order with cut of freq of 0.2
[b, a] = butter(6, 0.2);

%sampling time
time = 1/Fs: 1/Fs: numOfBits/dataRate;

%carrier
deltaf=.5;
f1 = Fc + (Fc*deltaf);
f2 = Fc - (Fc*deltaf);
carrier1 = cos(2 * pi * f1 * time);
carrier2 = cos(2 * pi * f2 * time);

%input
input = randi([0, 1], [1, 1024]);

%sampling rate is larger than data rate
%need to extend 1s and 0s by the ratio of the sampling rate and the data
%rate
extension_vector = ones(1, Fs/dataRate);
sampled_input = kron(input, extension_vector);
inverted_input = ~sampled_input
%for normal input multiply  by 1st carrier 
sampled_BFSK1 = sampled_input .* carrier1;
%inverted input multiplied by 2nd carrier
sampled_BFSK2 = inverted_input .* carrier2;
% (SUMMING MODULE) add the two signals to get the sampled BFSK
sampled_BFSK = sampled_BFSK1+sampled_BFSK2;

SNR = 5;

noise_variance = 1 / 10^(SNR/10);
noise_std = sqrt(noise_variance);
noise = noise_std .* randn(1, 1024 * Fs/dataRate);
received_signal = sampled_BFSK + noise;

%%=====================================demodulation%%=========================================================================
demodulated_signal1 = received_signal .* carrier1;
demodulated_signal2 = received_signal .* carrier2;
demodulated_signal = demodulated_signal1 + demodulated_signal2;

filtered_signal = filtfilt(b, a, demodulated_signal);
%% decoding signal (unfinished)
decoded_signal = zeros(1,1024);
for i = 1:1:1024
    interested_signal = filtered_signal(1 /2 * Fs/dataRate + (i - 1) * Fs/dataRate);
    if interested_signal > 0
        decoded_signal(i) = 1;
    else
        decoded_signal(i) = 0;
    end
end 
%% decoding signal (unfinished)
    
decoded_output = kron(decoded_signal, extension_vector);
ts1 = timeseries(sampled_input,time);
ts1.Name = 'Data waveform';
subplot(6, 1, 1);
plot(ts1);
xlim([0 0.01]);
ylim([-2 2]);

ts2 = timeseries(sampled_BFSK,time);
ts2.Name = 'Modulated BFSK Signal';
subplot(6, 1, 2);
plot(ts2);
xlim([0 0.01]);
ylim([-2 2]);

ts3 = timeseries(demodulated_signal,time);
ts3.Name = 'Received Signal';
subplot(6, 1, 3);
plot(ts3);
xlim([0 0.01]);
ylim([-4 4]);

ts4 = timeseries(filtered_signal,time);
ts4.Name = 'Demodulated signal';
subplot(6, 1, 4);
plot(ts4);
xlim([0 0.01]);
ylim([-4 4]);

ts5 = timeseries(decoded_output,time);
ts5.Name = 'Decoded signal';
subplot(6, 1, 5);
plot(ts5);
xlim([0 0.01]);
ylim([-2 2]);





