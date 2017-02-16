function [OPT] = problem_generator2(n,m,p_set,p_error,m1)
% This function is from Yuxin Chen http://web.stanford.edu/~yxchen/

if nargin < 5
    m1 = m;
end

stateDims = ones(1,n)*m;

ids = [1];
for i = 1:n
    ids1 = find(rand(1,m) > 1 - p_set);
    ids = [ids, ids1 + (i-1)*m+1];
    stateDims(i) = length(ids1);
end
h = 10;

OPT = linear_cons(m1,stateDims);

OPT.stateDims = stateDims;
OPT.W = [0, zeros(1, n*m);
    zeros(n*m,1), random_input(n,m, p_error)];
OPT.X_gt = [m, ones(1, n*m);
    ones(n*m,1), kron(ones(n,n), eye(m))];

OPT.W = OPT.W(ids, ids);
OPT.X_gt = OPT.X_gt(ids, ids);

function [W] = random_input(n, m, p_error)
dim = n*m;
W = eye(dim);

for i = 1:n
    for j = 1:n
        if i == j
            continue;
        end
        Wij = zeros(m, m);
        for k = 1:m
            if rand(1,1) > p_error
                Wij(k,k) = 1;
            else
                off = floor((m-1)*rand(1,1)) + 1;
                if off >= k
                    off = off + 1;
                end
                Wij(k, off) = 1;
            end
        end
        W(((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)) = Wij;
    end
end
W = (W+W')/2;

function [OPT] = linear_cons(m, stateDims)

dim = sum(stateDims)+1;
offsets = ones(1, length(stateDims)+1);

for i = 1:length(stateDims)
    offsets(i+1) = offsets(i) + stateDims(i);
end

% diagonal blocks
ids = 2:dim;
rowsA = 1:(2*dim-1);
colsA = [1, ((ids-1)*dim + ids), ((ids-1)*dim +1)];
valsA = ones(1, 2*dim-1);
b = [m, ones(1, 2*dim-2)];
numEqCons = 2*dim-1;

for i = 1:length(stateDims)
    for j = (offsets(i)+1):offsets(i+1)
        for k = (j+1):offsets(i+1)
            numEqCons = numEqCons + 1;
            id = (k-1)*dim + j;
            colsA = [colsA, id];
            rowsA = [rowsA, numEqCons];
            valsA = [valsA, 1];
            b = [b, 0];
        end
    end

end
OPT.A = sparse(rowsA, colsA, valsA, length(b), dim*dim);
OPT.b = b';

% diagonal blocks
rowsA = [1];
colsA = [1];
valsA = [1];
b = [m+1];
rowId = 1;
for i = 1:length(stateDims)
    for j = (i+1):length(stateDims)
        for k = (offsets(i)+1):offsets(i+1)
            rowId = rowId + 1;
            for l = (offsets(j)+1):offsets(j+1)
                id = (l-1)*dim + k;
                colsA = [colsA, id];
                rowsA = [rowsA, rowId];
                valsA = [valsA, 1];
            end
            b = [b, 1];
        end
    end
end
OPT.B1 = sparse(rowsA, colsA, valsA, length(b), dim*dim);
OPT.d1 = b';

rowsA = [];
colsA = [];
valsA = [];
b = [];
rowId = 0;
for i = 1:length(stateDims)
    for j = (i+1):length(stateDims)
        for l = (offsets(j)+1):offsets(j+1)
            rowId = rowId + 1;
            for k = (offsets(i)+1):offsets(i+1)
                id = (l-1)*dim + k;
                colsA = [colsA, id];
                rowsA = [rowsA, rowId];
                valsA = [valsA, 1];
            end
            b = [b, 1];
        end
    end
end
OPT.B2 = sparse(rowsA, colsA, valsA, length(b), dim*dim);
OPT.d2 = b';