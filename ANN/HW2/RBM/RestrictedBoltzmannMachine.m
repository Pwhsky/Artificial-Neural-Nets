
clear all
close all
clc

% Intialize XOR patterns
x = ...
    [-1 -1 -1;...
    1, -1, 1;...
    -1, 1, 1;...
    1, 1, -1;...
    1, 1, 1;...
    -1, -1, 1;...
    1, -1, -1;...
    -1, 1, -1];

% Initialize variables
N = 3;
M = 2; 
%[1 2 4 8];
eta = 0.01;


% Iteration variables
counter = 0;
nTrials = 100;
minibatchSize = 20;
k = 200;
N_out = 3000;
N_in = 2000;


V = zeros(N,1);
h = zeros(M,1);
theta_v = zeros(N,1);
theta_h = zeros(M,1);
w = normrnd(0,1, [M N]);


while true
    for trial = 1:nTrials
        dW = 0;
        dTheta_h = 0;
        dTheta_v = 0;

        for minibatch = 1:minibatchSize
            mu = randi(4);
            pattern = x(mu,:);
            v0 = pattern;
            b0_h = w*v0'-theta_h;
            Pb0_h = 1 ./ (1-exp(-2.*b0_h)); 

            for i = 1:M
                r = rand();
                if r < Pb0_h(i)
                    h(i) = 1;
                else
                    h(i) = -1;
                end
            end
            
            for t = 1:k %CD_k loop
                b_v = w'*h - theta_v;
                Pb_v = 1 ./ (1-exp(-2.*b_v));

                for j = 1:M
                    r = rand();
                    if r < Pb_v(j)
                        V(j) = 1;
                    else
                        V(j) = -1;
                    end
                end
                
                b_h = w*V - theta_h;
                Pb_h = 1 ./ (1-exp(-2.*b_h));

               
                for i = 1:M
                    r = rand();
                    if r < Pb_h(i)
                        h(j) = 1;
                    else
                        h(j) = -1;
                    end
                end
            end

           %Error computation and learning rule
           dW = dW + eta*(tanh(b0_h)*v0 - tanh(b_h)*V');
           dTheta_h = dTheta_h - eta*(v_0-V);
           dTheta_v = dTheta_v - eta*(tanh(b_h0) - tanh(b_h));

        end
        
        %Update weights and thresholds
        w = w + dW;
        theta_h = theta_h + dTheta_h;
        theta_v = theta_v + dTheta_v;

    end
    
    for i = 1:N_out
        indexPattern = randi(height(x));
        V = x(indexPattern,:);
        b_h = w*V - theta_h;
        
        Pb_h = 1 ./ (1-exp(-2.*b_h));

        for j = 1:M
            r = rand();
            if r < Pb_h(j)
                h(j) = 1;
            else
                h(j) = -1;
            end
        end

     
        for p = 1:N_in
        %%% here
        end
    end

end








