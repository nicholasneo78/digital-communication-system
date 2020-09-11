clear all; close all; clc;
%Define Signal length
%Assume the number of bits for transmission is 1024
Num_Bit = 1024;
%Define signal power                   
Signal_Power = 1; 
%SNR_dB = 10 log (Signal_Power/Noise_Power)                 
SNR_dB = 0:1:20;
%==> SNR = Signal_Power/Noise_Power = 10^(SNR_dB/10)
SNR = (10.^(SNR_dB/10));

%Set run times
Total_Run = 20;
%Different SNR value
for i = 1 : length(SNR)
 Avg_Error = 0;
 for j = 1 : Total_Run
  %Input singal
  %Generate random binary digits(0 or 1)
  Data = round(rand(1,Num_Bit));
  %Convert binary digit to (-1 or +1)
  Signal = 2 .* Data - 1;
  %Generate equal number of noise samples    
  Noise_Power = Signal_Power ./SNR(i);
        %The randn() function is for normal distribution
  Noise = sqrt(Noise_Power/2) .*randn(1,Num_Bit);   
  %Received Signal 
  Receive = Signal+Noise;                      

  Threshold = 0;
  Error = 0;
        %Fix the threshold value as 0
        %If received signal >= threshold value, threshold = 1
        %If received signal < threshold value, threshold = 0
  for k= 1 : Num_Bit
   if (Receive(k)>= Threshold) && Data(k)==0||(Receive(k)<Threshold && Data(k)==1)
    Error = Error+1;
   end
  end
  %Calculate bit error rate during transmission
  Error = Error ./Num_Bit;  
  %Calculate the average error for every runtime  
  Avg_Error = Error + Avg_Error;                   
 end
 Error_Rate(i) = Avg_Error / Total_Run;
end

%Calculate analytical Bit Error Rate
Theroy_Rate=(1/2)*erfc(sqrt(SNR)); 

%Graph and Plot the result           
figure(1)
semilogy (SNR_dB,Error_Rate,'b*');
ylabel('Pe');
xlabel('Eb/No')
hold on
semilogy (SNR_dB,Theroy_Rate,'k');
legend('simulation','theory',3);
axis([0 20 10^(-5) 1]);
hold off

%data generation
figure(2) 
subplot(311)
plot(Signal);
title('Data_Generated')
%noise generation
subplot(312) 
plot(Noise);
title('Noise_Generated')
%received data generation
subplot(313) 
plot(Receive);
title('Received_Data')