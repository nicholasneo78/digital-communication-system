clear all; close all; clc;

Fc = 10000; 
Fs = Fc * 16;
dataRate = 1000;
noOfBits = 1024;
amplitude = 1;

sampStart = 1/(2 * Fs);
sampInterval = 1/Fs;
timeTaken = noOfBits/dataRate;
time = sampStart: sampInterval: timeTaken;

carrier = amplitude .* cos(2 * pi * Fc * time);

input = randi([0, 1], [1, noOfBits]);
%Simulation for different SNR values
SNR = 0:5:50; %Loop from 0 to 50 (in multiples of 5)

[b, a] = butter(6, 0.2);

ratio_fs_dataRate = Fs/dataRate;
extension = ones(1, ratio_fs_dataRate);
sampled_input = kron(input, extension);

sampled_ook = sampled_input .* carrier;

SNRvalues = zeros(1,11); 
AverageOOKError = zeros(1,11);

S=1;
error = 0;
numOfRuns = 20;

for j = 1:numOfRuns
    bitErrorRateOutput = zeros(1,11);
    arrayIndex=1;
    noOfErrors = 0;
    for i = SNR
        SNRvalues(arrayIndex) = i;
        
        N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
        outputSignal = awgn(sampled_ook,i,N); %i=SNR, N=noise variance
        
        demod_signal = outputSignal .* (2 * carrier);
        
        filter_signal = filtfilt(b, a, demod_signal);
        
        decoded_signal = zeros(1,1024);
        
        for count=1:noOfBits
            produced_signal = filter_signal(1 /2 * Fs/dataRate + (count - 1) * Fs/dataRate);
            if (produced_signal > 0.5)
                decoded_signal(count) = 1;
            else
                decoded_signal(count) = 0;        
            end
          
        end
        
        SNRvalues(arrayIndex)=i; %Store current SNR value into array
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(arrayIndex)= bitErrorRate;
        AverageOOKError(arrayIndex) = AverageOOKError(arrayIndex) + bitErrorRateOutput(arrayIndex);
        
        arrayIndex = arrayIndex + 1;
    end
end

AverageOOKError = AverageOOKError ./ numOfRuns;
semilogy(SNRvalues, AverageOOKError, 'k-*');
ylim([10^(-5) 10^1]);
xlim([0 50]);
hold on

S=1;
SNRvalues = zeros(1,11);
bitErrorRateOutput = zeros(1,11);
arrayIndex=1;
AverageBPSKError = zeros(1,11);

sampled_input_bpsk = 2 * sampled_input - 1;
sampled_bpsk = sampled_input_bpsk .* carrier;

for runs = 1:numOfRuns
    bitErrorRateOutput = zeros(1,11);
    arrayIndex=1;
    for i = SNR
        SNRvalues(arrayIndex) = i;
        
        N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
        outputSignal = awgn(sampled_bpsk,i,N); %i=SNR, N=noise variance

        demod_signal = outputSignal .* (2 * carrier);
        
        filter_signal = filtfilt(b, a, demod_signal);
        
        decoded_signal = zeros(1,1024);
        
        for count=1:noOfBits
            produced_signal = filter_signal(1 /2 * Fs/dataRate + (count - 1) * Fs/dataRate);
            if produced_signal > 0
                decoded_signal(count) = 1;
            else
                decoded_signal(count) = 0;        
            end
        end
        
        SNRvalues(arrayIndex)=i; %Store current SNR value into array
        bitErrorRate = calculate_error_rate(decoded_signal, input);
        bitErrorRateOutput(arrayIndex)=bitErrorRate;
        AverageBPSKError(arrayIndex) = AverageBPSKError(arrayIndex) + bitErrorRateOutput(arrayIndex);
        arrayIndex = arrayIndex +1;
    end
end
AverageBPSKError = AverageBPSKError ./ numOfRuns;
semilogy(SNRvalues, AverageBPSKError, 'c-*');
axis([0 50 -1 1]);
ylabel('Log 10 Bit Error Rate') ;
hold on
title('Bit Error vs SNR - OOK and BPSK');
legend({'y = AverageOOK','y= AverageBPSK'},'Location','northeast')
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;

function bitErrorRate = calculate_error_rate(input, tempInput)
    error = 0;
    noOfBits= 1024;
    for i = 1:1:noOfBits
        if input(i) ~= tempInput(i)
            error = error + 1;
        end
    end
    bitErrorRate = error/noOfBits;
end
