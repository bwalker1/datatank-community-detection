%fn = '/netscr/walkeb6/tmp2079767123.mat';
%load(fn);

list = cell(nt,1);
addpath(genpath('/Users/ben/Dropbox/GenLouvain-2.1/'));
for i = 0:(nt-1)
    varname = ['Var_' int2str(i)];
    A = eval(varname);
    N = length(A);
    % convert distance to similarity value
    for j=1:N
        for k=(j+1):N
            if 1 || A(j,k) < threshold
                A(j,k)=exp(-A(j,k)/100);
            else
                A(j,k) = 0;
            end
        end
    end
    list{i+1} = (A + A');
end

N = length(list{1});
T = length(list);

[B,mm] = multiord(list,gamma,omega);
[S,Q] = iterated_genlouvain(B,10000,1,1,1);
Q = Q/mm;
Var=reshape(S,N,T);