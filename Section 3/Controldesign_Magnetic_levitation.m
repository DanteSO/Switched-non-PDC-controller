clear all
clc
%% System Parameters
m=0.05; % Kg
g=9.8; % m/s2
lambda=0.46; %H
mu=2; %m^-1
kk=0.001; %Ns/m
zref=0.04; %m
iref2=(2*m*g*(1+mu*zref)^2)/(lambda*mu); %A^2
%% Nonlinear functions (|x|<0.3)
f1min= (g*mu*(mu*0.3+2*mu*zref+2))/((1+mu*(zref+0.3))^2);
f1max= (g*mu*(mu*-0.3+2*mu*zref+2))/((1+mu*(zref-0.3))^2);
f2max= (-lambda*mu)/(2*m*(1+mu*(zref+0.3))^2);
f2min= (-lambda*mu)/(2*m*(1+mu*(zref-0.3))^2);
%% Local models
A{1}=[0 1; f1min -kk/m];
A{2}=[0 1; f1min -kk/m];
A{3}=[0 1; f1max -kk/m];
A{4}=[0 1; f1max -kk/m];
B{1}=[0; f2min];
B{2}=[0; f2max];
B{3}=[0; f2min];
B{4}=[0; f2max];

r=4; % number of rules
rho=100;
%}
%% INITIALIZING THE DECISION MATRICES
n=length(A{1});
mm=length(B{1}(1,:));
for i=1:r
    for j=1:r
    P{i}=sdpvar(n,n,'symmetric');
    K{i}{j}=sdpvar(mm,n);
    end
end
 %% DEFINITIONS OF SOME MATRICES

LMIs=[P{1}>=0];
    for i=1:r
        LMIs=[LMIs, P{i}>=0];
    end

sum=zeros(n);
for i=1:r
      sum= P{i}+sum;
end

% Extended matrix \xi_{ijk}

    for i=1:r
        for j=1:r
            for k=1:r
                 
YYYY{i}{j}{k}=  rho*(r*P{k}-sum)+A{i}*P{j}+P{j}*A{i}'+B{i}*K{k}{j}+K{k}{j}'*B{i}';
             
            end
        end
    end

%% DESIGN OF THE SWTICHED NON-PDC CONTROLLER 

 for k=1:r
        for i=1:r
            for j=1:r
                if i<=j
                
                    LMIs=[LMIs, 
    
                            YYYY{i}{j}{k}+YYYY{j}{i}{k}<=0      ];
    
 
         end
    end

      end
  end

%% Solver Configuration
opts=sdpsettings('solver','sdpt3','verbose',0,'warning',0);
sol = optimize(LMIs,[],opts);
[p,d]=checkset(LMIs);
checkset(LMIs)
clear sum
if sum(p < 0) == 0
    disp('XXXX  FEASIBLE XXXX')

    P1=value(P{1})
    P2=value(P{2})
    P3=value(P{3})
    P4=value(P{4})

    K11=value(K{1}{1})
    K12=value(K{1}{2})
    K13=value(K{1}{3})
    K14=value(K{1}{4})
    K21=value(K{2}{1})
    K22=value(K{2}{2})
    K23=value(K{2}{3})
    K24=value(K{2}{4})
    K31=value(K{3}{1})
    K32=value(K{3}{2})
    K33=value(K{3}{3})
    K34=value(K{3}{4})
    K41=value(K{4}{1})
    K42=value(K{4}{2})
    K43=value(K{4}{3})
    K44=value(K{4}{4})
else
   disp('XXXX INFEASIBLE XXXX')
end

%% Simulink model link

IC=[0.15 2.5]; % Initial conditions
t=0:0.0001:0.6; % Simulated Time
sim('Magnetic_levitation',t);

subplot(2,2,1)
plot(t,z_t)
subplot(2,2,2)
plot(t,x_t)
subplot(2,2,3)
plot(t,switchingfunction)
subplot(2,2,4)
plot(t,MFs)