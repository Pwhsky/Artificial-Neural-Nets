%Compute a time series prediction using a reservoir network layout
clc
clear
%Load data
tic
A = readmatrix("training-set.csv");
B = readmatrix("test-set-9.csv");

%Initialize 
kI = eye(500).*0.01;
M = 500;
N = 3;

time_steps = 500;
W_in = randn(500,N)*sqrt(0.002);
W_reservoir = randn(500)*sqrt(2/500);

%initial states of reservoir neurons:
r = zeros(500,1);
%R is the mega-matrix storing the r-states for all training patterns.
R = zeros(500,length(A));

%training
for o = 1:(length(A)-1)

x = A(:,o);
R(:,o) = r(:);

%Update rule
r = tanh(W_reservoir*r + W_in*x);

end %End of training

%Calculate the output matrix
W_out = A*R' * (R*R' + kI)^(-1);

%Feed test data
for o = 1:(length(B)-1)

x = B(:,o);
R(:,o) = r(:);
%Update rule

r = tanh(W_reservoir*r + W_in*x);

end
O = W_out*r;


%Predict the next 500 steps of the test data
for t = 1:time_steps
    
    r = tanh(W_reservoir * r + W_in * O);
    O = W_out*r;

    components(:,t) = O;
    
end

plot3(components(1,:),components(2,:),components(3,:))
y_components = components(2,:);
csvwrite("prediction.csv",y_components);

toc