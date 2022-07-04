addpath('../build/src','../build/mex');

N = 15;
xadj = [  1  3  6  9 12 14 17 21 25 29 32 34 37 40 43 45];
adjncy = [ 2  6  1  3  7  2  4  8  3  5  9  4 10  1  7 11  2  6  8 12  3 ...
    7  9 13  4  8 10 14  5  9 15  6 12  7 11 13  8 12 14  9 13 15 10 14];
partGraphRecursive = METIS_PartGraphRecursive(N,1,xadj,adjncy,[],[],[],2,[],[],[]);
partGraphKway = METIS_PartGraphKway(N,1,xadj,adjncy,[],[],[],3,[],[],[]);

n = 5;
A = 2*eye(n)-diag(ones(n-1,1),-1)-diag(ones(n-1,1),1);
A = kron(kron(A,eye(n)),eye(n))+kron(kron(eye(n),A),eye(n))+kron(kron(eye(n),eye(n)),A);

[perm,iperm] = ReorderSparseMat(A);
Aperm = A(perm,perm);

% NE = 5;
% NN = 7;
% eptr = [1 4 7 10 14 17];
% eind = [1 2 3 2 3 4 1 3 5 2 3 4 5 4 6 7];

NE = 4;
NN = 6;
eptr = [1 4 7 10 14];
eind = [1 2 3 2 3 4 1 3 5 2 3 4 5];

[epartMeshDual,npartMeshDual] = METIS_PartMeshDual(NE,NN,eptr,eind,[],[],1,2,[],[])
[epartMeshNodal,npartMeshNodal] = METIS_PartMeshNodal(NE,NN,eptr,eind,[],[],2,[],[]);
