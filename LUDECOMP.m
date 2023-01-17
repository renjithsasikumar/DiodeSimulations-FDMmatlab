function [Variable delta] = LUDECOMP(a,b,c,Variable,n_max,rHs)

%% LU decomposition 
 alpha(1) = b(1);
    for i=2:n_max
        beta(i)=a(i)/alpha(i-1);
        alpha(i)=b(i)-beta(i)*c(i-1);
    end
    
  %% 

    v(1) = rHs(1);
    for i = 2:n_max
        v(i) = rHs(i) - beta(i)*v(i-1);
    end
     
%%

    temp = v(n_max)/alpha(n_max);
    delta(n_max) = temp - Variable(n_max);
    Variable(n_max)=temp;
    for i = (n_max-1):-1:1       %delta%
        temp = (v(i)-c(i)*Variable(i+1))/alpha(i);
        delta(i) = temp - Variable(i);
        Variable(i) = temp;
    end
end