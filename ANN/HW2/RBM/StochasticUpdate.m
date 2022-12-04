function h_prime = StochasticUpdate(local_field)
 for i = 1:length(local_field)
           threshhold(i) = 1/(1+exp(-2*local_field(i)));
           r = rand(1);
           
           if r < threshhold(i)
               h_prime(i) = 1;
           else
           h_prime(i) = -1;
           end   
        end


end

