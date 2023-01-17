function f = BER(x)

%f=x/(exp(x)-1);

% absoX=abs(x);
% 
% if (absoX == 0)
%    f=1;
%    return ;
% end
if(x<=-23.0259)
    f=-x;
    return ;
end
if((x>-23.0259) && (x<-0.0198))
    f=x/(exp(x)-1);
    return ;
end
if((x>=-0.0198)&&(x<=0.0213))
    f=1 - ((x/2)*(1 - ((x/6)*(1 - x*x/60))));	
    return ;
end
if((x>0.0213)&&(x<23.0259))
    f=(x*exp(-x))/(1-exp(-x));
    return ;
end
if((x>=23.0259)&&(x<=26.2952))
    f=(x*exp(-x));
    return ;
end
if(x>=26.2952)
	%else
    f=0;
    return ;
end
end
