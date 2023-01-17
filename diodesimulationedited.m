%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%     1D Drift Diffusion Model fora pn Diodes       
%%%%%%%%%%     Equilibrium and Non Equilibrium Solver        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
      
clear all;
close all;
clc;

% Defining the Fundamental and Material Constants %

q     = 1.602E-19;        % C or [J/eV]
kb    = 1.38E-23;         % [J/K]
eps   = 1.05E-12;         % This includes the eps  = 11.7 for Si [F/cm]
T     = 300;              % [K]
ni    = 1.5E10;           % Intrinsic carrier concentration [1/cm^3]
Vt    = kb*T/q;           % [eV]
RNc   = 2.8E19;           % This is 2.8e20 in the FORTRAN file
TAUN0 = 0.1E-6;           % Electron SRH life time
TAUP0 = 0.1E-6;           % Hole SRH life time
mun0   = 1500;            % Electron Mobility in cm2/V-s
mup0   = 1000;            % Hole Mobility in cm2/V-s
                        
                        
dEc = Vt*log(RNc/ni);

% Define Doping Values %

Na = 1E17;             % [1/cm^3]
Nd = 1E16;             % [1/cm^3]

% Calculate relevant parameters for the simulation %

Vbi = Vt*log(Na*Nd/(ni*ni));
W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd));     % [cm]
Wn  = W*sqrt(Na/(Na+Nd));                    % [cm]
Wp  = W*sqrt(Nd/(Na+Nd));                    % [cm]
Wone = sqrt(2*eps*Vbi/(q*Na));               % [cm]
E_p = q*Nd*Wn/eps;                           % [V/cm]
Ldn = sqrt(eps*Vt/(q*Nd));
Ldp = sqrt(eps*Vt/(q*Na));
Ldi = sqrt(eps*Vt/(q*ni));

% Calculate relevant parameters in an input file %
% Write to a file
save input_params.txt  Na Nd Vbi W Wn Wp E_p Ldn Ldp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Material_Constants    %Define some material constants
      
% Setting the size of the simulation domain based 
% on the analytical results for the width of the depletion regions
% for a simple pn-diode %

x_max = 0;
if(x_max < Wn)
    x_max = Wn;
end
if(x_max < Wp)
    x_max = Wp;
end
x_max = 20*x_max;

% Setting the grid size based on the extrinsic Debye lengths %

dx = Ldn;
if(dx > Ldp)
    dx=Ldp;
end
dx = dx/20;

% Calculate the required number of grid points and renormalize dx %

n_max = x_max/dx;
n_max = round(n_max);

dx = dx/Ldi;    % Renormalize lengths with Ldi

% Set up the doping C(x)=Nd(x)-Na(x) that is normalized with ni %

for i = 1:n_max
    if(i <= n_max/2)
        dop(i) = - Na/ni;
    elseif(i > n_max/2)
        dop(i) = Nd/ni;
    end
end
% Initialize the potential based on the requirement of charge
% neutrality throughout the whole structure

for i = 1: n_max
    zz = 0.5*dop(i);
    if(zz > 0)
        xx = zz*(1 + sqrt(1+1/(zz*zz)));
    elseif(zz <  0)
        xx = zz*(1 - sqrt(1+1/(zz*zz)));
    end
    fi(i) = log(xx);
    n(i) = xx;
    p(i) = 1/xx;
end

delta_acc = 1E-5;               % Preset the Tolerance


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%           EQUILIBRIUM  SOLUTION PART BEGINS                      %%

%(A) Define the elements of the coefficient matrix for the internal nodes and
%    initialize the forcing function

dx2 = dx*dx;
for i = 1: n_max
    a(i) = 1/dx2;
    c(i) = 1/dx2;
    b(i) = -(2/dx2+exp(fi(i))+exp(-fi(i)));
    f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i))+exp(-fi(i)));
 end


%(B) Define the elements of the coefficient matrix and initialize the forcing
%    function at the ohmic contacts 

a(1) = 0;
c(1) = 0;
b(1) = 1;
f(1) = fi(1);
a(n_max) = 0;
c(n_max) = 0;
b(n_max) = 1;
f(n_max) = fi(n_max);

%(C)  Start the iterative procedure for the solution of the linearized Poisson
%     equation using LU decomposition method:

flag_conv = 0;		           % convergence of the Poisson loop
k_iter= 0;
while(~flag_conv)            
    k_iter = k_iter + 1; 
    
    alpha(1) = b(1);
    for i=2:n_max
        beta(i)=a(i)/alpha(i-1);
        alpha(i)=b(i)-beta(i)*c(i-1);
    end
    
% Solution of Lv = f %    

    v(1) = f(1);
    for i = 2:n_max
        v(i) = f(i) - beta(i)*v(i-1);
    end
     
% Solution of U*fi = v %    

    temp = v(n_max)/alpha(n_max);
    delta(n_max) = temp - fi(n_max);
    fi(n_max)=temp;
    for i = (n_max-1):-1:1       %delta%
        temp = (v(i)-c(i)*fi(i+1))/alpha(i);
        delta(i) = temp - fi(i);
        fi(i) = temp;
    end
    
    delta_max = 0;
    
    for i = 1: n_max
        xx = abs(delta(i));
        if(xx > delta_max)
            delta_max=xx;
        end
        %sprintf('delta_max = %d',delta_max)      %'k_iter = %d',k_iter,'
        
    end

     %delta_max=max(abs(delta));
   
% Test convergence and recalculate forcing function and 
% central coefficient b if necessary
    
    if(delta_max < delta_acc)
        flag_conv = 1;
    else
        for i = 2: n_max-1
            b(i) = -(2/dx2 + exp(fi(i)) + exp(-fi(i)));
            f(i) = exp(fi(i)) - exp(-fi(i)) - dop(i) - fi(i)*(exp(fi(i)) + exp(-fi(i)));
        end
    end
end


xx1(1) = dx*1e4;
for i = 2:n_max-1 
    Ec(i) = dEc - Vt*fi(i);     %Values from the second Node%
    ro(i) = -ni*(exp(fi(i)) - exp(-fi(i)) - dop(i));
    el_field1(i) = -(fi(i+1) - fi(i))*Vt/(dx*Ldi);
    el_field2(i) = -(fi(i+1) - fi(i-1))*Vt/(2*dx*Ldi);
    n(i) = exp(fi(i));
    p(i) = exp(-fi(i));
    xx1(i) = xx1(i-1) + dx*Ldi*1e4;
end


Ec(1) = Ec(2);
Ec(n_max) = Ec(n_max-1);
xx1(n_max) = xx1(n_max-1) + dx*Ldi*1e4;
el_field1(1) = el_field1(2);
el_field2(1) = el_field2(2);
el_field1(n_max) = el_field1(n_max-1);
el_field2(n_max) = el_field2(n_max-1);
nf = n*ni;
pf = p*ni;
ro(1) = ro(2);
ro(n_max) = ro(n_max-1);

figure(1)
plot(xx1, Vt*fi,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Potential [eV]');
title('Potential vs Position - at Equilibrium');

figure(2)
plot(xx1, el_field1,'r','LineWidth',2)
hold on;
plot(xx1, el_field2,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Electric Field [V/cm]');
title('Field Profile vs Position - at Equilibrium');

figure(3)
%plot(xx1, nf,'g','LineWidth',2)
semilogy(xx1, nf,'g','LineWidth',2)
hold on;
%plot(xx1, pf,'r','LineWidth',2)
semilogy(xx1, pf,'r','LineWidth',2)
xlabel('x [um]');
ylabel('Electron & Hole Densities [1/cm^3]');
title('Electron & Hole Densities vs Position - at Equilibrium');
legend('n','p');
%axis([0 6.75 0 10.2e17])

figure(4)
%plot(xx1, ro,'r','LineWidth',2)
plot(xx1, q*ro,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Total Charge Density [C/cm^3]');
title('Total Charge Density vs Position - at Equilibrium');
%axis([0.5 5 -3e17 8e17])

figure(5)
plot(xx1, Ec,'r','LineWidth',2)
xlabel('x [um]');
%ylabel('Total Charge Density [1/cm^3]');
ylabel('Conduction Band Energy (eV)');
title('Conduction Band vs Position - at Equilibrium');