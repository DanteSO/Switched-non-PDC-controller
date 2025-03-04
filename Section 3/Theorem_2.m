clear all
clc
%% STEPS TO OBTAIN THE FEASIBILITY REGION 

    for b=0:0.25:8
    for a=30:1:60

%% EXAMPLE_1

r=4; % number of rules
rho=2;

A{1}=[1 -4;-1 -2];
A{2}=[-3 -a;1 -2];
A{3}=[1 -4;20 -2];
A{4}=[-3 -4;20 -2];

B{1}=[1;20];
B{2}=[1;20];
B{3}=[1;-b/10];
B{4}=[1;1];

%% INITIALIZING THE DECISION MATRICES
n=length(A{1});
m=length(B{1}(1,:));
for i=1:r
    for j=1:r
    P{i}=sdpvar(n,n,'symmetric');
K{i}{j}=sdpvar(m,n);
R1{i}= sdpvar(n,n,'full');
R2{i}= sdpvar(n,n,'full');
    end
end

%% DEFINITIONS OF SOME MATRICES

% Parametric matrices P_{i}
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
                 
 YYYY{i}{j}{k}=  [rho*(r*P{k}-sum)+A{i}*P{j}+P{j}*A{i}'+B{i}*K{j}{k}+K{j}{k}'*B{i}'];
               
            end
        end
    end

%% SET THE LMIs OF THEOREM 2

 for k=1:r
        for i=1:r
            for j=1:r
                if i<=j
                
                    LMIs=[LMIs, 
    
                            YYYY{i}{j}{k}+YYYY{j}{i}{k}<=0;      ];
    
 
         end
    end

      end
  end

%% SOLVER CONFIGURATION
opts=sdpsettings('solver','sdpt3','verbose',0,'warning',0);
sol = optimize(LMIs,[],opts);
[p,d]=checkset(LMIs);
checkset(LMIs)
clear sum
if sum(p < 0) == 0
    disp('XXXX  FEASIBLE XXXX')

% Put a circle if the LMI is feasible    
plot(a,b,'o','MarkerFaceColor','white','MarkerEdgeColor','k')
hold on

else
   disp('XXXX INFEASIBLE XXXX')
end
 end 
end 
