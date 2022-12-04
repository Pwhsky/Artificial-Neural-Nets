%Restricted Boltzmann Machine
%3 visible, 1,2,4,8 hidden neurons with signum function.

%k is the iteration ammount

clear
data_set = [[-1,-1,-1]',[1,-1,1]',[-1,1,1]',[1,1,-1]'];
M_set = [1,2,4,8];
epochs = 200;
iterations = 1e2;
learnRate = 0.03;




for o = 1:length(M_set)
    M = M_set(o);

    
%Initialize weights & theta 
    
W  = 0.2*randn(M,3);
theta_v = 0.2*(randn(3,1));
theta_h = 0.2*(randn(M,1));




for e = 1: epochs
        
%We have 4 patterns to choose from
    
mu = data_set(:,randi(4));
v_j = mu;

%initialize the visible neurons to the pattern:
%initialize hidden neurons
v_h = W*v_j - theta_h;

for i = 1:M
     r = rand;
    P = (1+exp(-2*v_h(i)))^-1;
    
    if r<P
        v_h(i) = 1;
    else
        v_h(i) = -1;
    end
    
end
v_0j = v_j;
v_0h = v_h;




%update visible neurons:
v_j = W'*v_h -  theta_v;

for i = 1:3
    r = rand;
    P = (1+exp(-2*v_j(i)))^-1;
    
    if r<P
        v_j(i) = 1;
    else
        v_j(i) = -1;
    end
    
end


%Update hidden neurons:
v_h = W*v_j - theta_h;

for i = 1:M
     r = rand;
    P = (1+exp(-2*v_h(i)))^-1;

    if v_h(i)<r
        v_h(i) = 1;
    else
        v_h(i) = -1;
    end
    
end

%Update rule:
%In essence, update using the difference between the initiated value and
%the last value. 

dw = learnRate*(tanh(v_0h)*v_0j' - tanh(v_h)*v_j');
dtheta_v = learnRate*(v_0j - v_j);
dtheta_h = learnRate*(tanh(v_0h) - tanh(v_h));

W = W+dw;
theta_v = theta_v +dtheta_v;
theta_h = theta_h + dtheta_h;


end

e = exp(0.5*  W* v_j*v_j')


                              
end
%Calculate Kullback-Leibler divergence:


%KL_d(o) =  sum(0.25.*log(0.25./ ));


