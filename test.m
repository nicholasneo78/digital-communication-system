% simple matlab script from notes
% semi colon at the back prevents printing of variable values
% press f5 to run script
a=5;
b=a/2;
k=3*a;
a
b
k
k^2

% all variables are matrices in MATLAB
% row vector
row = [12,14,-1];
row
% column vector
col = [13;42;-5];
col

% matrix (, and ;)
matrix = [1,2,3;4,5,6;7,8,9];
matrix

% sub matrix
% sub_matrix = matrix(r1:r2,c1:c2)
col_2nd = matrix(:,2);
col_2nd
row_2nd = matrix(2,:);
row_2nd

% simple if-else statement
k = rand*100; % rand is a single uniformly distributed random number in the interval (0,1)
% check value
k
if k <=30
    "0 to 30 range"
elseif k<=60
    "30 to 60 range"
else
    "60 to 100 range"
end

% for loop
number = 0
for i = 1:10
   number = number+i;
end

% while loop
number2 = 0
i = 0;
while i<10
    number2 = number2 + i;
    i = i+1;
end
number
number2

% Scalar Matrices
a = 3;
b = [1,2,3;4,5,6];
c = b+a;
d = b-a;
e = b*a;
f = b/a;
c
d
e
f

% change the numbers if you want to see how the curves will change
t = 0:0.01:1;
f = 5;
y = cos(2*pi*f*t);
plot(t,y)

title("Cosine Wave");
xlabel("Time in seconds");
ylabel("Amplitude");

% plot sine and cosine
t = 0:0.01:1;
f = 5;
y1 = sin(2*pi*f*t);
y2 = cos(2*pi*f*t);
subplot(2,1,1)
plot(t,y1);
title("Sine Wave");
xlabel("Time in seconds");
ylabel("Amplitude");
subplot(2,1,2)
plot(t,y2);
title("Cosine Wave");
xlabel("Time in seconds");
ylabel("Amplitude");