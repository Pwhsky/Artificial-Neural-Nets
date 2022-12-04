%Homework 1:
%Calculation of error probabilities for asynchronous updates
%in a Hopfield Network.     

%Author: Alex Lech
%Date of completion: 2022-09-05


clc
clear
number_of_bits = 120;

patterns_vector = [12,24,48,70,100,120];

iterations = 1000;

%///Occurence_vector records number of errors & output_probabilities is the
%calculated chance of error.
occurence_vector = zeros(6,1);
output_probabilities = zeros(1,6);
%--------------------------------------------------------------------------

Hebbs = false;
Hebbs_type = input("Set diagonal weights to 0?      Y/N:","s");
if Hebbs_type == "Y"
   Hebbs = true ;
else
    Hebbs = false;
end
tic


for o = 1:length(patterns_vector)
    storage_size = patterns_vector(o) ;
   % pattern_storage_matrix = zeros(number_of_bits,storage_size);
    errors = 0;
    
for y = 1:iterations
 

%///---Random patterns creation
pattern_storage_matrix =  randi([0 1],number_of_bits,storage_size);
pattern_storage_matrix(pattern_storage_matrix==0)=-1;

%for j = 1:storage_size
     %for i = 1:number_of_bits 
     %   if pattern_storage_matrix(i,j) == 0
     % pattern_storage_matrix(i,j) = -1;
     %   end
     %end
%end


%zero_list = any(pattern_storage_matrix,'all');
%pattern_storage_matrix(zero_list) = -1;




%///----------------------------


%///Generate weight matrix W for a random pattern v

v = randi([1 storage_size]);
W =  zeros(120,120);
for mu = 1:storage_size 
       if mu ~= v   
       Weights=( pattern_storage_matrix(1:120,mu)*pattern_storage_matrix(1:120,mu)');
       W = W+Weights; 
       end

end   
W = W*1/120;

%///-----------------------------------------------


%///Code block to set diagonal weights to 0, this is used in task 1
if Hebbs
    
        for i = 1:number_of_bits

             W(i,i) = 0;
           
        end
end
%-------------------------------------------------------------------
    
    
        

  %///Neuron update, using equation (2.28)
  random_neuron = randi([1 120]);
  activation =  sign((1-1/120)*pattern_storage_matrix(random_neuron,v) + W(random_neuron, 1:120)*pattern_storage_matrix(1:120,v));
      
  if activation == 0
  activation = 1;
  end
  %///-------------------------------------
  
  
  %///did the random_neuron change?
  if activation ~= pattern_storage_matrix(random_neuron,v)
      errors = errors + 1;
  end
  %---------------------------------
  
end
  
%///record the number of errors and calculate error probability.
occurence_vector(o) = errors;
output_probabilities(o) = errors/(iterations);
output_probabilities

%----------------------------------------------------------------

end
total_time = toc