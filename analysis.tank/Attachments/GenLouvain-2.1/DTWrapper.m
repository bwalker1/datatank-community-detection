sz = size(Var);
Var = reshape(Var,sqrt(sz(1)),sqrt(sz(1)),sz(2));

list = cell(nt,1);
addpath(genpath('/nas02/w/a/walkeb6/GenLouvain-2.1/'));
for i = 1:nt
    A = squeeze(Var(:,:,i));
    N = length(A);
    % convert distance to similarity value
    for j=1:N
        for k=1:N
            if (j==k)
                continue;
            end

            if A(j,k) < threshold
                A(j,k)=exp(-A(j,k)/100);
            else
                A(j,k) = 0;
            end
        end
    end
    list{i} = (A);
end

N = length(list{1});
T = length(list);

[B,mm] = multiord(list,gamma,omega);
[S,Q] = iterated_genlouvain(B,10000,1,1,1);
Q = Q/mm;
Var=reshape(S,N,T);
