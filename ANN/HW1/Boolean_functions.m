
clc
clear
tic
dimension = [2,3,4,5];

learnRate = 0.05;
iterations = 10000;
epochs= 20;
n = 2;

%for n = 2:2  


 if n>3
      max_length = iterations;
    else
        max_length = 2^(2^n);
 end
 
 error_list = zeros(1,max_length);
 
 function_list = strings(1,max_length);
 
 

 


for y = 1:max_length
    
%Sample a target function.    
target =  randi([0 1],1,2^n);
target(target==0)=-1;

%Create a string copy to append to a list.
copy_string = strjoin(string(target));

    
 %Check if it's in the list, resample untill a unique pattern is found 
    while ismember(copy_string,function_list) == true &&  length(function_list) < max_length
    target =  randi([0 1],1,2^n);
    target(target==0)=-1;
    copy_string = strjoin(string(target));
    function_list(y) = copy_string;
    end
    

    function_list(y) = copy_string;
    %intiialize a weight vector
W = (n^(-0.5))*randn(1,n);

%Initialize targets.
inputs = (dec2bin(0:(2^n)-1)' - '0');
inputs(inputs==0) =-1;


%Initialize Threshhold

theta = zeros(1,2^n);

%theta = 0;

    
  

dW = 0;
for e = 1:epochs
       
%Output term:

for i = 1:2^n
O(i) =  W*inputs(:,i) - theta(i);
end


O(O >=0)=1;
O(O <0)=-1; 


%Update weights and theta:
for i = 2^n
  x_j= inputs(:,i);
  
dW = dW + learnRate.*(target(i)-O(i))*x_j;

theta(i) = theta(i) -learnRate.*(target(i) - O(i));
end

W = W+ dW';

 
%Break loop if all outputs are mapped to target.
correct = 0;
for i = 1:2^n
if O(i) == target(i)
    correct = correct +1;
end
end

if correct == 2^n
    break
end



end %End of epoch loop

error = 0;

%Check if O = t
for i = 1:2^n
if O(i) ~= target(i)
    error = error +1;
    
end
end


%One can check the index for an error occurence in error_list,
%and then locate it in function_list. 

error_list(y) = error

%Append the unique function to the list.



end

%n_of_errors(n) = sum(error_list)


%end
toc
