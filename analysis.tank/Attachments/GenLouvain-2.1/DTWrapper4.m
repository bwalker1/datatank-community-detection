load('Input.mat');
addpath(genpath('HelperFunctions/'));
Var = adjacency;
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

[B,mm] = multiord(list,gamma,omega);
tic
[S,Q] = genlouvain(B,10000,1,1,1);
toc
Q = Q/mm;
Var=reshape(S,N,T);

Var = Var;
save('Output.mat','Var','-v4');
