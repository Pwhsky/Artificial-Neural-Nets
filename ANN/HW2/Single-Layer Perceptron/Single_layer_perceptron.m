%Single-Layer Perceptron
%Task: Experiment with different number of neurons in hidden layer
%minimize error-classification to below 12%.

clear
tic
A = readmatrix('training_set');
B = readmatrix('validation_set');

learnRate = 0.03;


% Set the activation function and its derivative
g = @(b) tanh(b);
gprim = @(b) (1-(tanh(b).^2));


tMax = 1e6; % Number of time steps
epochs = 15;
% Store patterns in P_training, set mean to 0 and variance to 1

P_train = A(:,[1:2])';
meanP = mean(P_train,2);
varP = sqrt(var(P_train')');
P_train = (P_train - meanP)./varP;
B(:,1:2) = (B(:,1:2) - meanP')./varP';
targets = A(:,3)';




%% Train the perceptron:

% Initialize & pre-allocate
O = zeros(1,length(P_train));


    K = 12;
    
    % Initialize weights
    w = 0.1*(randn(K,2)-1);
    W = 0.1*(randn(1,K)-1);
    theta = zeros(K,1)-1;
    Theta = zeros(1,1)-1;
    
    
    for t = 0:tMax-1
        

        % Choose a random pattern
        mu=randi(length(A));
        xi = P_train(:,mu);
        
        for e = 1:epochs
        
        % Forward-propagation
        b_j = w*xi - theta;
        V = g(b_j);
        
        b_i = W*V - Theta;
        O(mu) = g(b_i);
        
        % Error
        d_i = gprim(b_i).*(targets(mu)-O(mu));
        
        d_j = (d_i'*W').*gprim(b_j);
        
        % Back-propagation
      %  w = w + learnRate*d_j*xi';
      
        w = w + learnRate*d_j*xi';
        
        W = W + learnRate*d_i*V';
        
        theta = theta - learnRate*d_j;
        Theta = Theta - learnRate*d_i ;
        end
   end
    %Classification error:
    C = 1/(2*length(P_train))*sum(abs(sign(O)-targets));
     
    %Energy function:
    H = 1/2 * sum( (targets - O).^2 );
    
%To get the desired format of the weights, transpose.
writematrix(w,'w1.csv')
writematrix(W','w2.csv')

writematrix(theta,'t1.csv')
writematrix(Theta,'t2.csv')



%41 attempts on openTA has to be some kind of record...

toc
