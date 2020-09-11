%[1-1]
%Generate random binary digits (0 or 1) 
noOfBits = 1024;
inputData = randi([0, 1], [1, noOfBits]);

%[1-2]
%Convert binary digits to +/- 1 (1->+1, 0->-1) => Data for transmission
%Replace 0 to -1
inputData(inputData==0) = -1;

%[1-14]
%Simulation for different SNR values
SNR = 0:5:50; %Loop from 0 to 50 (in multiples of 5)

%[1-6]
%Assume input data has unit power
S=1; 

%[1-16]
bitErrorRate = zeros(1,11); %Initialise an array of 1-by-11 zeros to store error rate
SNRvalues = zeros(1,11); %Initialise an array of 1-by-11 zeros to store SNR values
arrayIndex = 1; %To count the array index for storage

for i = SNR %Loop from 0 to 50 (in multiples of 5)
    %[1-11]
    noOfErrors = 0; %Count number of errors for each SNR
    
    %[1-6]
    N=S./(10.^(i./10)); %Obtain noise variance (10log10 = S/N)
    
    %[1-3 && 1-7]
    %AWGN function applys White Gaussian Noise to inputData
    outputSignal = awgn(inputData,i,N); %i=SNR, N=noise variance
    
    %[1-8]
    %Convert outputSignal to 0/1
    %(outputSignal >= 0) -> 1
    %(outputSignal < 0) -> 0
    threshold = 0; 
    outputSignal(outputSignal>=threshold) = 1;
    outputSignal(outputSignal<threshold) = 0;
    
    %[1-9]
    %Conpute bit error rate
    for count=1:noOfBits %Loop from 1 to 1024
        %If inputData > 0 but outputSignal < 0, error
        %If inputData <= 0 but outputSignal >= 0, error
        if (inputData(count) > threshold && outputSignal(count) == 0) || (inputData(count) <= threshold && outputSignal(count) == 1)
            noOfErrors = noOfErrors + 1;
        end
    end
    %[1-16]
    SNRvalues(arrayIndex)=i; %Store current SNR value into array
    bitErrorRate(arrayIndex)= noOfErrors/noOfBits; %Store the bitErrorRate into array
    arrayIndex=arrayIndex+1;
end

%[1-15]
%Plot the results
semilogy(SNRvalues, bitErrorRate);
axis([0 50 -1 1])
title("Plot of Bit Error Rate vs Signal to Noise Ratio");
xlabel('E_{b}/N_{0}') ;
ylabel('P_{e}') ;
