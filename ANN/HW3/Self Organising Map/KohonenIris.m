clear
tic
A = readmatrix("iris-data.csv");
B = readmatrix("iris-labels.csv");

sigma_0 = 10;
decay_sigma = 0.05;

eta_0 = 0.1;
decay_eta = 0.01;

%standardise input data:
input_data = A./(max(max(A)));

W = rand([40,40,4]);
W_original = W;
epochs = 10;

output = zeros(40);

%Neighborhood function
h = @(distance,sigma) exp( -((distance).^2) ./(2*sigma.^2));


for e = 1:epochs
           
    sigma = sigma_0*exp(-decay_sigma * e);
    eta = eta_0*exp(-decay_eta * e);
    
    three_sigma = 3*sigma;
   for p = 1:length(input_data)
       %sample a pattern:
       rIndex = randi(length(input_data));
       
       x = input_data(rIndex,:);
       xlabel = B(rIndex);
       
       %Sample a random weight
      
       
       %calculate distance between weight and pattern
       
       distance = W()
        for i = 1:40
            for j = 1:40
                
                r = [i j];
                distance_r0 = vecnorm(r-r0);

             if distance_r0 < three_sigma     
             W(i,j,:) = W(i,j,:) +  eta*h(distance_r0,sigma)*(x'-W(i,j,:));
             end
    end
  end
   
   end
end
toc

function [i0 j0] = WinningNeuron(x,W)



end

hold on