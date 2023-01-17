%%%%%%%%%%%%  Code has been done for Numericalmethode applied to
%%%%%%%%%%%%  Problem1(DIFF equations)
%%%%%%%%%%%%  It uses forward and backward difference---can include
%%%%%%%%%%%%  central difference using conditional statements while
%%%%%%%%%%%%  creating w matrix and rHs..................... %%%%%%%%%%%%






close all;
clear all;
clc;

%%%%%%%%% Reading the input parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 P = input('Enter the value of applied force');
 l = input('Enter the length of the bar');
 A = input('Enter the area of cross section');
 q = input('Enter the value of q');
 x=0:.1:l;
 E = 2*10^5;
 
 %%%%%%%%%%% Calculating displacement and Sigma using Exact Equation %%%%%%%%%%%
 
 displacement = (P/(E*A))*x + (q/(E*A))*(x.*(l-(x/2)));
 sigma= (P/A)+(q/A)*(l-x);
 
 %pause;
%%%%%%%%%%%%%%      Creat Matrix   %%%%%%%%%%%%%%%
 
 a=[1 0 0 ];
 b=[0 -1 1];
 c=[1 -2 1];
w=[]; 
n = input('Enter the number of nodes');
w(1,:)=[a zeros(1,n-3)];
w(2,:)=[c zeros(1,n-3)];
w(n,:)=[zeros(1,n-3) b];

for i= 3:n-1
    h= w(i,:);
    h(i-1:i+1)=c;
    w(i,:)=h;
end
%%%%%%%%%%%%%%      Creat rHs Matrix   %%%%%%%%%%%%%%%
atForward = (-q*(l^2))/((n-1)^2*E*A);
atBackward = ((P*l)/((n-1)*E*A));

rHs=[];
rHs(1,:)=0;
rHs(n,:)=atBackward;
rHs(2:n-1)=atForward;

%%%%%%%%%%%%% Find and plot  Displacement using  Numerical methode %%%%%
displacementNumerical = inv(w)*rHs;
y=0:(l/(n-1)):l
plot(x,displacement,'r--',y,displacementNumerical,'y--');
title('Exact Vs FDM ................backward ');
xlabel('length');
ylabel('Displacement');
legend('Exact','FDM');

 for i=1:(n-1)
     sigmaNumerical(i)= E*((displacementNumerical(i+1)-displacementNumerical(i))/(l/(n-1)));
     z(i) = (y(i+1)+y(i))/2;
 end
 figure;
 plot(x,sigma,'r--',z,sigmaNumerical,'b-');
 title('Exact Vs FDM ................backward ');
xlabel('length');
ylabel('Stress');
legend('Exact','FDM');

 
 
     
     
     
     
     