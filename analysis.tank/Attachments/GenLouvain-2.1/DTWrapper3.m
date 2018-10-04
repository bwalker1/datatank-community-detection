load('Input.mat');
addpath(genpath('HelperFunctions/'));
Var = network;
sz = size(Var);
Var = reshape(Var,sqrt(sz(1)),sqrt(sz(1)),sz(2));
network = Var;
%quit
nt = sz(2);
list = cell(nt,1);
for i = 1:nt
    A = squeeze(Var(:,:,i));
    list{i} = (A);
end

N = length(list{1});
T = length(list);

coeffs = [];
%save('Output.mat','N','-v4');
for run=1:1
    if runs > 1
        gamma = gmin + (gmax-gmin)*rand;
        omega = wmin + (wmax-wmin)*rand;
    else
        gamma = gmin;
        omega = wmin;cd 
    end
    [B,mm] = multiord(list,gamma,omega);
    tic
    [S,Q] = genlouvain(B,10000,1,1,1);
    toc
    Q = Q/mm;
    Var=reshape(S,N,T);
    
    
    % compute the coefficients for champ
    A=0;P=0;C=0;
    %disp(size(network));
    edgeWeightSum = squeeze(sum(network));
    sliceWeightSum = squeeze(sum(edgeWeightSum));
    sz = size(network);
    for t=1:T
        for i=1:N
            for j=1:N
                if i==j
                    continue;
                end
                if Var(i,t)==Var(j,t)
                    A = A + network(i,j,t);
                    P = P + edgeWeightSum(i,t)*edgeWeightSum(j,t)/sliceWeightSum(t);
                end
            end
            if t < T-1
                if Var(i,t)==Var(i,t+1)
                    C = C + 1;
                end
            end
        end
    end
    coeffs = [coeffs;A,P,C];
end
Var = Var;
save('Output.mat','Var','-v4');
