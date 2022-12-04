clear
tic
train_data = [[-1 -1 -1]',[1 -1 1]', [-1 1 1]', [1 1 -1]'];

full_dataset = [[-1 -1 -1]',[1 -1 1]', [-1 1 1]', [1 1 -1]',[1,1,1]',[-1,-1,1]', [1,-1,-1]',[-1,1,-1]'];



N = 3;
hidden_neurons = [1,2,4,8];


%Parameters: 
learnRate = 0.005;

batches = 1;
k= 1000;
trials = 1000;
N_out = 3000;
N_inner = 3000;
normalized_probability = (N_out*N_inner)^(-1);

DKL = zeros(length(hidden_neurons),1);
DKL_est = zeros(length(hidden_neurons),1);

P_B=zeros(length(hidden_neurons),8);


pattern_frequency_matrix= zeros(length(hidden_neurons),length(full_dataset));

mb = 40;
P_D = 1/4;



for i = 1:length(hidden_neurons)
    
    M = hidden_neurons(i);
    

%Initialize weights & threshholds
W = normrnd(0,1/sqrt(N),[M,N]);
theta_h = zeros(M,1);
theta_v = zeros(1,N);

%initialize neurons:
v = zeros(N,1); 
h = zeros(M,1);


%Training
for t = 1:trials
%Initialize weights & threshholds

dW = zeros(M,N);
dtheta_v = zeros(1,N);
dtheta_h  = zeros(M,1);  





    %v_0 and h_0 will be local fields, whereas h and v are -1 or 1.
    
for m = 1:mb
    
 pattern_index = randi(4);
 pattern = train_data(:,pattern_index);
 
 %initialize the visible neurons value to pattern
   v_0 = v; 
   v = pattern;
  
    %Update hidden neurons, start by evaluating the field
    h_field = (W*v_0-theta_h); 
    h_0 = h_field;
    %Stochastic update of h:
    h = StochasticUpdate(h_field);
       
  
    for j = 1:k %Cdk loop
    %Update all visible neurons: compute local field and perform a stochastic update.
    
    v_field = h*W -theta_v;
    v = (StochasticUpdate(v_field))';
    
    %Update hidden neurons, like above:
    h_field = W*v-theta_h;
    h = StochasticUpdate(h_field);
    
    end
    
    
    dW    = dW +  learnRate*(tanh((h_0))*(v_0') - tanh(h_field)*v') ;
    dtheta_v = dtheta_v +  learnRate*( (v_0' - v') );
    dtheta_h = dtheta_h + learnRate*(tanh(h_0) - tanh(h_field));
end
W = W +dW;
theta_v = theta_v + dtheta_v;
theta_h = theta_h + dtheta_h;


end



%training ends here!

  %%
            
%select random pattern in the full list of functions, will it be recognised?

for j = 1:N_out
   
    random_index = randi(8);
    random_pattern = full_dataset(:,random_index);
    v = random_pattern;
    
    %evaluate the field
    b_h = W*v-theta_h; 
    h = StochasticUpdate(b_h);
  
    for y = 1:N_inner
        %evaluate the field of visible neurons:
        b_v = h*W-theta_v;
        v = StochasticUpdate(b_v)';
        
        
        
         %evaluate the field of hidden neurons:
         b_h = W*v-theta_h;
         h = StochasticUpdate(b_h);
        
  
        %record the number of recognised patterns:
     
        if v == random_pattern
           P_B(i,random_index) = P_B(i,random_index) + normalized_probability;
       end
        
  

    end

end


%for u = 1:4
 
        
DKL_est(i) =   DKL_est(i) + P_D *  log(P_D/(sum(P_B(i,1:4))) );


    
%end


  
   
end



  
    
hold on

plot(hidden_neurons,DKL_est)

    
%%

 % Dkl theoretical limit.
DKL(1) = (3 - log2(1) - (1)/(2^(log2(1))));
for i = 2:8
   
    
if i < 2^(N-1)-1
    DKL(i) = 3 - log2(i+1) - (i+1)/(2^(log2(i+1)));
else
    DKL(i) =0;
    
end
end
hold on
plot([1:8],log(2)*DKL)
toc
xlabel("Neurons")
ylabel("KL-Divergence")
legend("Calculated DKL","Theoretical limit")


function output = StochasticUpdate(local_field)
 for i = 1:length(local_field)
           threshhold(i) = 1/(1+exp(-2*local_field(i)));
           r = rand(1);
           
           if r < threshhold(i)
               output(i) = 1;
           else
           output(i) = -1;
           end   
        end


end



