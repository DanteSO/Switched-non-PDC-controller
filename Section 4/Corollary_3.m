clear all
clc
 %% EXAMPLE 5

a=20; % change the values to see the feasibility
b=5.30; % change the values to see the feasibility

r=2; % number of rules
rho=1;

A{1}=[3.6 -1.6;6.2 -4.3];
A{2}=[-a -1.6;6.2 -4.3];
polyA = rolmipvar(A,'A',r,1);
B{1}=[-0.45;-3];
B{2}=[-b;-3];
polyB = rolmipvar(B,'B',r,1);

%% INITIALIZING THE DECISION MATRICES

n = length(A{1});
N=4; %Polinomal degree

indices = cell(1, N);  % N is the number of indices
[indices{:}] = ndgrid(1:r);
indices = cell2mat(cellfun(@(x) x(:), indices, 'UniformOutput', false));

%% INITIALIZING THE DECISION MATRICES

for row = 1:size(indices, 1)
    idx = num2cell(indices(row, :));
    for k = 1:r
        P{idx{:}} = sdpvar(n, n, 'symmetric');
        K{k}{idx{:}} = sdpvar(1, n);
    end
end

%% REDUCING THE COMPUTATIONAL COMPLEXITY BASED ON REMARK 9

f=0;
for i=1:r
    cont{i}=0;
end


for row = 1:size(indices, 1)
    idx = num2cell(indices(row, :));
    v = cell2mat(cellfun(@(x) x(:), idx, 'UniformOutput', false));
             f=f+1;
             PP{f}={v_ind(v,r),P{idx{:}}};
        for k=1:r
            KK{k}{f}={v_ind(v,r),K{k}{idx{:}}};
            if sign(n_rep(v,k))==1
                vv{k}=v;
                cont{k}=cont{k}+1;
                g=find(v == k);
                vv{k}(g(1))=[];
% We eliminate the redundante matrices P{i} substituing by P_{i}
                P_{k}{cont{k}}={v_ind(vv{k},r),n_rep(v,k)*P{idx{:}}};
            end
        end            
 end

 
polyPP = rolmipvar(PP,'PP',r,N);

for k=1:r
    polyKK{k} = rolmipvar(KK{k},'KKk',r,N);
    polyP{k} = rolmipvar(P_{k},'PPP',r,N-1);
end


%% DEFINITIONS OF SOME MATRICES

LMIs=[polyPP>=0];

sum1=zeros(n);
for k=1:r
      sum1= polyP{k}+sum1;
end

% Extended matrix \xi_{ijk}

    for k=1:r
                 
        YYYY{k}=  [rho*(r*polyP{k}-sum1)+polyA*polyPP+polyPP*polyA'+...
                   polyB*polyKK{k}+polyKK{k}'*polyB'];
                
    end


%% SET THE LMIs OF COROLLARY 3

 for k=1:r
 
         LMIs=[LMIs,  YYYY{k}<=0];

 end

%% SOLVER CONFIGURATION
opts=sdpsettings('solver','sdpt3','verbose',0,'warning',0);
sol = optimize(LMIs,[],opts);
[p,d]=checkset(LMIs);
checkset(LMIs)
clear sum
if sum(p < 0) == 0
    disp('XXXX  FEASIBLE XXXX')
    P1=value(P{1});
    P2=value(P{2});
else
   disp('XXXX INFEASIBLE XXXX')
end





