%% Drift- Diffussion solver using Gummels method Version 2.1 %%%%%%%%%%%%%%%%%%%%%%%%
%%% using FDM - as explained in Prof. D.Vasileska ASU %%%%%%%%%%%%%%%%%%%
%%% Version: 2. 0 Equillirium solver alone%      %%%%%%%%%%%%%%%%%%%%%%%%
%%% Version2.1 added the Non-Equillibrum solver %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;
clear all;
close all;
clc;
%%
% Defining the Fundamental and Material Constants %

q     = 1.602E-19;        % C or [J/eV]
kb    = 1.38E-23;         % [J/K]
eps0  = 8.854E-14;        % [F/cm]
epsr  = 11.7;             % Enter the relative permitivity of the material
eps   = eps0*epsr;        % [F/cm]
T     = 300;              % [K]
ni    = 1.5E10;           % Intrinsic carrier concentration [/cm^3]
Vt    = kb*T/q;           % [eV]
RNc   = 2.8E19;           % *****
TAUN0 = .01E-6;           % Electron SRH life time
TAUP0 = .01E-6;           % Hole SRH life time
mun0   = 1410;            % Electron Mobility [cm2/V-s]
mup0   = 470;            % Hole Mobility in [cm2/V-s]
                        
                        
dEc = Vt*log(RNc/ni);

% Define Doping Values %

Na = 1E16;             % [/cm^3]
Nd = 1E17;             % [/cm^3]

% Calculate relevant parameters for the simulation %

Vbi = Vt*log(Na*Nd/(ni*ni));                % Built in potential[V]
W   = sqrt(2*eps*(Na+Nd)*Vbi/(q*Na*Nd));    % Depletion width [cm]
Wn  = W*sqrt(Na/(Na+Nd));                   % Depletion width in n side[cm]
Wp  = W*sqrt(Nd/(Na+Nd));                   % Depletion width in P side[cm]
Wone = sqrt(2*eps*Vbi/(q*Na));              % [cm]
E_p = q*Nd*Wn/eps;                          % Peak Electric field[V/cm]
Ldn = sqrt(eps*Vt/(q*Nd));                  % Debay lenth [cm]
Ldp = sqrt(eps*Vt/(q*Na));                  % Debay lenth [cm]
Ldi = sqrt(eps*Vt/(q*ni));                  % intrinsic Debay lenth [cm]

%%

% Define the device dimensions

       devDim  = 400E-4;                     % Total device deminention in [cm]
       x_fracP=0.5;                         % Deines the length of Pside [0-0.9]
%%
% Meshing settings :  Meshing Starts here
       %dx = min(Ldn, Ldp);
       %dx=  1.5000e-05;
       dx=devDim/200;
      % n_max=round(devDim/dx);
       n_max = 200;
       dx = dx/Ldi;                         % Renormalize lengths with Ldi
%% 
% Define the doping profile: Here we define the device structre  
       dopCon(1:ceil(x_fracP*n_max))= -Na/ni;
       dopCon(ceil(x_fracP*n_max)+1:n_max)= Nd/ni;
%% n=0.5*(C+math.sqrt(C*C+4*ni*ni))
% Ensuring the charge nuetrality on entire device: 
        for i = 1: n_max
        zz = 0.5*dopCon(i);
        if(zz > 0)                           % N side (1 give result than 4ac)
            xx = zz*(1 + sqrt(1+1/(zz*zz)));
        elseif(zz <  0)                      % P side
            xx = zz*(1 - sqrt(1+1/(zz*zz)));
        end
% Gives the initial potential profile from carrier concentation eeuation
% Express n,p as function of Potential, fermi potential,intrisic.concentation and Ei, Fermi
% potential[two equations in note: compare: corrected equations in box]
% Ei(considerinfg as the refrence levl).
% Assume 
        Phi_old(i) = log(xx);                     
        n(i) = xx;                      
        p(i) = 1/xx;
        end

        delta_acc = 1E-5;                    % Preset the Tolerance

%%  %%%%%%%%%%%%%%%%%%% Equillibrium solver begin here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Initially assume all node based on doping profile.
% 2. Calculate the coeffcient matrix.
        % 2.1 Discretize using FDM
        % 2.2 Phi_new,Coeffecient matrix,Phi_old
        % 2.3 a,b,c are three elements of Sparse matrx row.
        % 2.4 Discretization steps in note and Vasileska FDM slide
% 3. Solve AX=B. for Phi_new
% 4. Do convergencs check [ ] Cov check commented while
% loop
% 5. Check for tolerance. 
% 6. If converge stop : or ugo to step 3.
%%
    dx2 = dx*dx;                             % Calculating dxsquare     
%   Define the coefficents for the matrix 
    flag=0;
    jk=0;
    
    
    while( ~ flag)
        jk=jk+1
    for i = 2: n_max-1
        a(i) = 1/dx2;
        c(i) = 1/dx2;
        b(i) = (-2/dx2-(exp(Phi_old(i))+exp(-Phi_old(i)))); 
        rHs(i) = exp(Phi_old(i)) - exp(-Phi_old(i)) - dopCon(i) - Phi_old(i)*(exp(Phi_old(i))+exp(-Phi_old(i)));
    end
%   Deine the boundary conditions.
    a(1) = 0;
    c(1) = 0;
    b(1) = 1;
    rHs(1) = Phi_old(1);
    rHs(n_max) = Phi_old(n_max);
    a(n_max) = 0;
    c(n_max) = 0;
    b(n_max) = 1;
%   Boundary conditions end here.

%%
    % Solve the matrix using LU decomposition.
    [Phi_old delta]=LUDECOMP(a,b,c,Phi_old,n_max,rHs);
    deltaPhi_max=(max(abs(delta)));
    Phi_matrix_old(:,jk)= Phi_old;
 
%% Check convergence

    if(deltaPhi_max< delta_acc )
        flag =1;       
    end
%%    
    end;
%%%%%%%%%%%%%%%%% Equillibrium solver ends here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% Plot the results. %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%     1. Do the enrgy level shift
%     2. Calculate lectric field profils, carrier density
%     3. Do the denormalization steps.
    
xx1(1) = dx*1e4*Ldi;
for i = 2:n_max-1 
    Ec(i) = dEc - Vt*Phi_old(i);     %Values from the second Node%
    ro(i) = -ni*(exp(Phi_old(i)) - exp(-Phi_old(i)) - dopCon(i));
    el_field1(i) = -(Phi_old(i+1) - Phi_old(i))*Vt/(dx*Ldi);
    el_field2(i) = -(Phi_old(i+1) - Phi_old(i-1))*Vt/(2*dx*Ldi);
    n(i) = exp(Phi_old(i));
    p(i) = exp(-Phi_old(i));
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
% 
equi_Phi_old = Phi_old;
% Calculate Quasi Fermi Level - Efn Efp
        for i = 1:n_max
            Ei(i)   = Ec(i) - 0.56;
            Efn(i)  = Ei(i) + Vt*log(nf(i)/ni);
            Efp(i)  = Ei(i) - Vt*log(pf(i)/ni);
        end
        Ev = Ec - 1.12;                         % material property Eg=1.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot the results. %%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% plot(xx1, Vt*Phi_old,'r','LineWidth',2)
% xlabel('x [um]');
% ylabel('Potential [eV]');
% title('Potential vs Position - at Equilibrium');
% 
% figure(2)
% plot(xx1, el_field1,'r','LineWidth',2)
% hold on;
% plot(xx1, el_field2,'r','LineWidth',2)
% xlabel('x [um]');
% ylabel('Electric Field [V/cm]');
% title('Field Profile vs Position - at Equilibrium');
% 
% figure(3)
% plot(xx1, nf,'g','LineWidth',2)
% semilogy(xx1, nf,'g','LineWidth',2)
% hold on;
% plot(xx1, pf,'r','LineWidth',2)
% semilogy(xx1, pf,'r','LineWidth',2)
% xlabel('x [um]');
% ylabel('Electron & Hole Densities [1/cm^3]');
% title('Electron & Hole Densities vs Position - at Equilibrium');
% legend('n','p');
% %axis([0 6.75 0 10.2e17])
% 
% figure(4)
% %plot(xx1, ro,'r','LineWidth',2)
% plot(xx1, q*ro,'r','LineWidth',2)
% xlabel('x [um]');
% %ylabel('Total Charge Density [1/cm^3]');
% ylabel('Total Charge Density [C/cm^3]');
% title('Total Charge Density vs Position - at Equilibrium');
% %axis([0.5 5 -3e17 8e17])
% 
% figure(5)
% plot(xx1, Ec,'r','LineWidth',2)
% xlabel('x [um]');
% %ylabel('Total Charge Density [1/cm^3]');
% ylabel('Conduction Band Energy (eV)');
% title('Conduction Band vs Position - at Equilibrium');

%% Plot Energy band diagrams with quasi fermi levels.
        figure(11)
        plot(xx1, Ec,'black','LineWidth',2.5);
        hold on;
        plot(xx1, Ev,'black','LineWidth',2.5);
        hold on;
        plot(xx1, Ei,'--g','LineWidth',2.5);
        hold on;
        plot(xx1, Efn,'r','LineWidth',2.5);
        hold on;
        plot(xx1, Efp,'b','LineWidth',2.5);
        xlabel('x [um]');
        ylabel('Energy [eV]');
        title('Quasi Fermi Levels (Efn & Efp) vs Position - at Applied Bias(0.625V)');
        legend('Ec','Ev','Ei','Efn','Efp');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);


        
% %% Non equillibrium solution starts here: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Define the mobility profile of the device  %%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Assuming constant mobility: Neglecting, Doping, Teperature and High field
% % effects on mobility.
vindex=0;

%% Normalize mu amd Tau %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mubar = max(mun0,mup0);
mun0 = mun0/mubar;
mup0 = mup0/mubar;
TAUBAR = (Ldi*Ldi)/(mubar*Vt);
TAUN0=TAUN0/TAUBAR;
TAUP0=TAUP0/TAUBAR;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% setting  convergence parametes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delta_acc = 1E-5;		% COVEGRGENCE MARGIN
alpha=0.005; 			% DAMPING PARAMETER
%Phi_old_int = Phi_old(1);
%% Volatge step increments.
for VA = -0.2:.33*Vt:.625
    flag2=0;
    vindex = vindex +1
    if(vindex==38)
        check=1;
    end
    
    Vplot(vindex) = VA;
    Phi_old(1) = Phi_old(1)+ VA; 
    %% Initialize the First and Last Node for Poisson's eqn
    a(1) = 0;
    c(1) = 0;
    b(1) = 1;
    a(n_max) = 0;
    c(n_max) = 0;
    b(n_max) = 1;
    rHs(1) = Phi_old(1);
    rHs(n_max) = Phi_old(n_max); 
    for i = 1:n_max           
        mup(i) = mup0;               
        mun(i) = mun0;                    
    end 
    
while(~flag2)
%%
%% Asuume constant mobility at every node.
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% 3.2 Solve Continuity Equation for Electron and Holes                                    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Boundary conditions @ Anode
        fp(1) = p(1);
        fn(1) = n(1);
        an(1) = 0;              
        bn(1) = 1;              
        cn(1) = 0;              
        ap(1) = 0;              
        bp(1) = 1;              
        cp(1) = 0;              
        % Boundary conditions @ Cathode
        fn(n_max) = n(n_max);
        fp(n_max) = p(n_max);
        an(n_max) = 0;          
        bn(n_max) = 1;          
        cn(n_max) = 0;          
        ap(n_max) = 0;          
        bp(n_max) = 1;          
        cp(n_max) = 0;    
      
   for i = 2: n_max-1
            
            munim1by2 = (mun(i-1)+mun(i))/2;         
            munip1by2 = (mun(i)+mun(i+1))/2;         
            mupim1by2 = (mup(i-1)+mup(i))/2;         
            mupip1by2 = (mup(i)+mup(i+1))/2;
            %% Co-efficients for HOLE Continuity eqn
            cp(i) = mupip1by2 * BER(Phi_old(i) - Phi_old(i+1));
            ap(i) = mupim1by2 * BER(Phi_old(i) - Phi_old(i-1));
            bp(i) = -( mupim1by2 * BER(Phi_old(i-1) - Phi_old(i)) + mupip1by2 * BER(Phi_old(i+1) - Phi_old(i)));
	    %% Forcing Function for HOLE Continuity eqns
            fp(i) = (Ldi*Ldi*dx2/(Vt*TAUBAR*mubar)) * ( p(i)*n(i) - 1 ) / ( TAUP0*(n(i) + 1) + TAUN0*(p(i)+1)); % SRH
            %fp(i) = (Ldi*Ldi*dx2/(Vt*TAUBAR*mubar)) * ( p(i) - 1 ) /TAUP0;  % Direct
        end
      %% Solve continuity equations p.  
            pold=p;
            [p deltap] = LUDECOMP(ap,bp,cp,p,n_max,fp); 
            p= alpha*p+(1-alpha)*pold;                 % DAMPING
   for i = 2: n_max-1
            
            %% Co-efficients for ELECTRON Continuity eqn
            cn(i) = munip1by2 * BER(Phi_old(i+1) - Phi_old(i));
            an(i) = munim1by2 * BER(Phi_old(i-1) - Phi_old(i));
            bn(i) = -( munim1by2 * BER(Phi_old(i) - Phi_old(i-1)) + munip1by2 * BER(Phi_old(i) - Phi_old(i+1)));

            %% Forcing Function for ELECTRON  Continuity eqns
            fn(i) = (Ldi*Ldi*dx2/(Vt*TAUBAR*mubar)) * ( p(i)*n(i) - 1 ) / ( TAUP0*(n(i) + 1) + TAUN0*(p(i)+1)); % SRH
                     
              %fn(i) = (Ldi*Ldi*dx2/(Vt*TAUBAR*mubar)) * ( n(i) - 1 ) /TAUN0;  % Direct
              %fp(i) = (Ldi*Ldi*dx2/(Vt*TAUBAR*mubar)) * ( p(i) - 1 ) /TAUP0;  % Direct
            
        end    
%% Solve continuity equations n. 
 nold=n;                 
[n deltan] = LUDECOMP(an,bn,cn,n,n_max,fn); 
n= alpha*n+(1-alpha)*nold;
%% Tracking n and p
       nmatrix(:,vindex)=n;
       pmatrix(:,vindex)=p;
%%   
% % Recalculate forcing function and central coefficient b for fi 
% % here values of n(i) and p(i) are used in place of exp(fi(i))
    
        for i = 2: n_max-1
              b(i) = -(2/dx2 + n(i) + p(i));
              rHs(i) = n(i) - p(i) - dopCon(i) - (Phi_old(i)*(n(i) + p(i)));
              %b(i) = -(2/dx2);
              %rHs(i) = n(i) - p(i) + dopCon(i) ; 
        end
%% Solve Poissons equation.
    Phi_oldn = Phi_old;
    [Phi_old delta]=LUDECOMP(a,b,c,Phi_old,n_max,rHs);
    deltaPhi_max=(max(abs(delta)));
    Phi_matrix(:,vindex)= Phi_old;
   
Phi_old= alpha*Phi_old+(1-alpha)*Phi_oldn;
%% Convergence  check
        if(deltaPhi_max< delta_acc )
         flag2 =1;
        end;


end
for i=2:n_max-1
        Jnim1by2(vindex,i) = (q*mun(i)*Vt/(dx*Ldi)) * ni*( n(i)*BER(Phi_old(i)-Phi_old(i-1)) - n(i-1)*BER(Phi_old(i-1)-Phi_old(i)) );
        Jnip1by2(vindex,i) = (q*mun(i)*Vt/(dx*Ldi)) * ni*( n(i+1)*BER(Phi_old(i+1)-Phi_old(i)) - n(i)*BER(Phi_old(i)-Phi_old(i+1)) );
        Jelec(vindex,i) = (Jnip1by2(vindex,i) + Jnim1by2(vindex,i))/2;
        Jpim1by2(vindex,i) = (q*mup(i)*Vt/(dx*Ldi)) * ni*( p(i)*BER((Phi_old(i-1)-Phi_old(i))) - p(i-1)*BER((Phi_old(i)-Phi_old(i-1))) );
        Jpip1by2(vindex,i) = (q*mup(i)*Vt/(dx*Ldi)) * ni*( p(i+1)*BER((Phi_old(i)-Phi_old(i+1))) - p(i)*BER((Phi_old(i+1)-Phi_old(i))) );
        Jhole(vindex,i) = (Jpip1by2(vindex,i) + Jpim1by2(vindex,i))/2;
 end

        Jtotal = (Jelec + Jhole)*mubar;

end
% Since the cuurrent is constant througth the device 
        Jtotal(:,1) = Jtotal(:,2);
        Jelec(:,1) = Jelec(:,2);
        Jhole(:,1) = Jhole(:,2);
        Jtotal(:,n_max) = Jtotal(:,(n_max-1));
        Jelec(:,n_max) = Jelec(:,(n_max-1));
        Jhole(:,n_max) = Jhole(:,(n_max-1));
        
%%
            xx1(1) = dx*1e4*Ldi;
        for i = 2:n_max-1 
            Ec(i) = dEc - Vt*Phi_old(i);     
            ro(i) = -ni*(n(i) - p(i) - dopCon(i));
            el_field1(i) = -(Phi_old(i+1) - Phi_old(i))*Vt/(dx*Ldi);
            el_field2(i) = -(Phi_old(i+1) - Phi_old(i-1))*Vt/(2*dx*Ldi);
            xx1(i) = xx1(i-1) + dx*Ldi*1e4;
        end

%%  Calculate field and carrier densities at non equillibrium    
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

%% Calculate Quasi Fermi Level - Efn Efp
        for i = 1:n_max
            Ei(i)   = Ec(i) - 0.56;
            Efn(i)  = Ei(i) + Vt*log(nf(i)/ni);
            Efp(i)  = Ei(i) - Vt*log(pf(i)/ni);
        end
        Ev = Ec - 1.12;                         % material property Eg=1.2
%% Plot carrier distribution at V=0.625V
        figure(6)
        semilogy(xx1,nf,'g',xx1,pf,'r','LineWidth',2);
       %semilogy(xx1(3:n_max),(n(3:n_max)),'g',xx1(3:n_max),p(3:n_max),'r','LineWidth',2);
        xlabel('x [um]');
        ylabel('Carrier density [1/cm^3]');
        %ylabel('Conduction Band Energy (eV)');
        title('Elec. & Hole concentrations - at Applied Bias(0.625V)');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);
%% Plot the potential.
        figure(7)
        plot(xx1, Vt*Phi_old,'r',xx1, Vt*equi_Phi_old ,'b','LineWidth',2)
        xlabel('x [um]');
        ylabel('Potential [eV]');
        title('Potential vs Position - at Applied Bias(0.625V)');
        legend('@.625V','@Equillibrium');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);
%% Plot current
   
        figure(8)
        plot(Vplot, Jtotal(:,n_max),'r','LineWidth',2)
        xlabel('VA [V]');
        ylabel('Total Current Density [Amp/cm^2]');
        title('I vs V Plot');
        %legend('Jtotal','Jhole','Jelec','2');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);
%% Plot Energy band diagrams with quasi fermi levels.
        figure(9)
        plot(xx1, Ec,'black','LineWidth',2.5);
        hold on;
        plot(xx1, Ev,'black','LineWidth',2.5);
        hold on;
        plot(xx1, Ei,'--g','LineWidth',2.5);
        hold on;
        plot(xx1, Efn,'r','LineWidth',2.5);
        hold on;
        plot(xx1, Efp,'b','LineWidth',2.5);
        xlabel('x [um]');
        ylabel('Energy [eV]');
        title('Quasi Fermi Levels (Efn & Efp) vs Position - at Applied Bias(0.625V)');
        legend('Ec','Ev','Ei','Efn','Efp');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);
% %% Plot current
%    
        figure(10)
        semilogy(Vplot, Jtotal(:,n_max),'r','LineWidth',2)
        xlabel('VA [V]');
        ylabel('Total Current Density [Amp/cm^2]');
        title('IV (log)');
        legend('Jtotal','Jhole','Jelec','2');
        set(findobj(gca,'type','axes'),'FontSize',18,'FontWeight','Normal','LineWidth',0.6);
