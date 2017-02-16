clear all;
close all;
filelist = dir('.\');

% all noise in the mesh was added in this code error induced preventing it
% from running remove error k = 
L = 1;
for i = 4:numel(filelist)
a = filelist(i).name;
b = '.\';
b = strcat(b,a);
%p =struct([]);
if(strcmp(b(end-3:end),'.off'))
    fid=fopen(b);
    fgetl(fid);
    nos = fscanf(fid, '%d %d  %d', [3 1]);
    nopts = nos(1);
    notrg = nos(2);
    coord = fscanf(fid, '%g %g  %g', [3 nopts]);
    p(L).field = a;
    p(L).coord = coord';
    triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
    triang=triang';
    p(L).triang=triang(:,2:4)+1;
    
    L = L+1;
    fclose(fid);
end
%%we have added 1 because the vertex indices start from 0 in vtk format


end
%load('centaur_HKS.mat');
%figure(1);
%hold on
%plot3(p(1).coord(:,1),p(1).coord(:,2),p(1).coord(:,3),'g.');

%plot3(p(1).coord(1:200:3400,1),p(1).coord(1:200:3400,2),p(1).coord(1:200:3400,3),'r.');
%hold off
points = [1:125:12500];
load('hkslist');
s1 = P(1).hks;
s2 = P(2).hks;
s3 = P(3).hks;
s4 = P(4).hks;
s1 = s1(points,:);
s2 = s2(points,:);
s3 = s3(points,:);
s4 = s4(points,:);
Srel = zeros(400,400);
S1 = zeros(100,100,4);
S2 = zeros(100,100,4);
S3 = zeros(100,100,4);
S4 = zeros(100,100,4);

    for i = 1:100
       S1(:,i,1)  = sum((s1 - repmat(s1(i,:),100,1)).^2,2);
       S1(:,i,2)  = sum((s1 - repmat(s2(i,:),100,1)).^2,2);
       S1(:,i,3)  = sum((s1 - repmat(s3(i,:),100,1)).^2,2);
       S1(:,i,4)  = sum((s1 - repmat(s4(i,:),100,1)).^2,2);
    end
    for i = 1:100
       S2(:,i,1)  = sum((s2 - repmat(s1(i,:),100,1)).^2,2);
       S2(:,i,2)  = sum((s2 - repmat(s2(i,:),100,1)).^2,2);
       S2(:,i,3)  = sum((s2 - repmat(s3(i,:),100,1)).^2,2);
       S2(:,i,4)  = sum((s2 - repmat(s4(i,:),100,1)).^2,2);
    end
    for i = 1:100
       S3(:,i,1)  = sum((s3 - repmat(s1(i,:),100,1)).^2,2);
       S3(:,i,2)  = sum((s3 - repmat(s2(i,:),100,1)).^2,2);
       S3(:,i,3)  = sum((s3 - repmat(s3(i,:),100,1)).^2,2);
       S3(:,i,4)  = sum((s3 - repmat(s4(i,:),100,1)).^2,2);
    end
    for i = 1:100
       S4(:,i,1)  = sum((s4 - repmat(s1(i,:),100,1)).^2,2);
       S4(:,i,2)  = sum((s4 - repmat(s2(i,:),100,1)).^2,2);
       S4(:,i,3)  = sum((s4 - repmat(s3(i,:),100,1)).^2,2);
       S4(:,i,4)  = sum((s4 - repmat(s4(i,:),100,1)).^2,2);
    end
    dimGroup = [100 100 100 100];
    
Srel = [S1(:,:,1) S1(:,:,2) S1(:,:,3) S1(:,:,4);S2(:,:,1) S2(:,:,2) S2(:,:,3) S2(:,:,4);S3(:,:,1) S3(:,:,2) S3(:,:,3) S3(:,:,4);S4(:,:,1) S4(:,:,2) S4(:,:,3) S4(:,:,4)];
alpha = 0.2;
Wrel = -Srel + alpha*(ones(400,400));
W = sparse(Wrel);
[X,info,A] = mmatch_CVX_ALS(W,dimGroup);
X = full(X);
colorvec = zeros(1,100);
for i = 1:100
    colorvec(i) =  200*(i-1);
end
X11 = X(1:100,1:100);
X12 = X(1:100,101:200);
X13 = X(1:100,201:300);
X14 = X(1:100,301:400);
X21 = X(101:200,1:100);
X22 = X(101:200,101:200);
X23 = X(101:200,201:300);
X24 = X(101:200,301:400);
X31 = X(201:300,1:100);
X32 = X(201:300,101:200);
X33 = X(201:300,201:300);
X34 = X(201:300,301:400);
X41 = X(301:400,1:100);
X42 = X(301:400,101:200);
X43 = X(301:400,201:300);
X44 = X(301:400,301:400);


load('deodesic');
%%
addpath gm-toolbox
[N,M] = find(X23 == 1);
match = [N';M'];
n = length(N);
sim = zeros(1,n);
for i = 1:n
    sim(i) = S2(N(i),M(i),3);
end

wEdge = 10;
if wEdge <= 0
    Xraw = double(sim);
else
    % construct graphs
    fprintf('- constructing graphs ...\n');
    [uniq_feat1,~,new_feat1] = unique(match(1,:));
    [uniq_feat2,~,new_feat2] = unique(match(2,:));
    cand_matchlist_uniq = [new_feat1';new_feat2'];
    xy1 = (p(2).coord)';
    xy2 = (p(3).coord)';
    edgeAttr1 = computeEdgeAttr(xy1(:,uniq_feat1));
    edgeAttr2 = computeEdgeAttr(xy2(:,uniq_feat2));
    fE1 = makeFeatBin(edgeAttr1);
    fE2 = makeFeatBin(edgeAttr2);
    eSimVal = computeDotProdSimilarity_sym( int32(cand_matchlist_uniq), fE1, fE2); % symmetric affinities
    vSimVal = sim;
    mSimVal = repmat(vSimVal,[numel(vSimVal),1]); % distribute the v score
    affinityMatrix = wEdge*eSimVal + mSimVal + mSimVal';
    affinityMatrix(1:size(affinityMatrix,1)+1:end) = 0;
    affinityMatrix(affinityMatrix<1e-3) = 0;
    affinityMatrix = sparse(double(affinityMatrix));
    fprintf('- running random walk matching ...\n');
    [group1,group2] = make_group12(match(1:2,:));
    Xraw = RRWM(affinityMatrix,group1,group2);
end
X12distort = greedyMatch(match,Xraw);

k = find(X12distort);
match1 = match(:,k);

x23 = zeros(100,100);
[n,m] = size(match1);
for i = 1:m
x23(match1(1,i),match1(2,i)) = 1;
end

%%

v1 = p(1).coord;    t1 = p(1).triang;
v2 = p(3).coord;    t2 = p(2).triang;

n = size(v1, 1);

[X1_i, X2_i] = find(x13);
X1_i = points(X1_i);
X2_i = points(X2_i);

colors = 1:10:10*size(X1_i,2);

colorMap_X1 = zeros(n,1);
colorMap_X2 = zeros(n,1);

% assign he colors to the matched points given by mesh_index
colorMap_X1(X1_i) = colors;
colorMap_X2(X2_i) = colors;

for i=1:n
    if(colorMap_X1(i) == 0)
        [m, index] = min(Gdist1(i,X1_i));
        colorMap_X1(i) = colorMap_X1(X1_i(index));
    end
    if(colorMap_X2(i) == 0)
        [m, index] = min(Gdist3(i,X2_i));
        colorMap_X2(i) = colorMap_X2(X2_i(index));
    end
end
%%

v3 = p(3).coord;    t3 = p(3).triang;
v4 = p(4).coord;    t4 = p(4).triang;
n1 = size(v3, 1);
[X3_i, X4_i] = find(x34);
X3_i = points(X3_i);
X4_i = points(X4_i);
colors1 = sort(randperm(1024, size(X3_i,2)));

colorMap_X3 = zeros(n1,1);
colorMap_X4 = zeros(n1,1);

colorMap_X3(X3_i) = colors1;
colorMap_X4(X4_i) = colors1;

for i=1:n1
    if(colorMap_X3(i) == 0)
        [m, index] = min(Gdist3(i,X3_i));
        colorMap_X3(i) = colorMap_X3(X3_i(index));
    end
    if(colorMap_X4(i) == 0)
        [m, index] = min(Gdist4(i,X4_i));
        colorMap_X4(i) = colorMap_X4(X4_i(index));
    end
end

%%
figure(1)
scatter3(v1(:,1), v1(:,2), v1(:,3), 3, colorMap_X1)
title('mesh1 ColorMap');

figure(2)
scatter3(v2(:,1), v2(:,2), v2(:,3), 3, colorMap_X2)
title('mesh2 ColorMap');

%%
figure(3)
scatter3(v3(:,1), v3(:,2), v3(:,3), 3, colorMap_X3)
title('mesh3 ColorMap');

figure(4)
scatter3(v4(:,1), v4(:,2), v4(:,3), 3, colorMap_X4)
title('mesh4 ColorMap');


%% plot corresponding lines

cm = jet(size(X1_i,2));
idx = floor(1+ (colors - min(colors)) / (max(colors) -min(colors)) * (size(X1_i,2) - 1) );
rgb_mat = cm(idx,:);

v2_shift = v2;
v2_shift(:,2) = v2(:,2) + (max(v2(:,2)) - min(v2(:,2))+0.2);

figure(2)
hold on
scatter3(v1(:,1), v1(:,2), v1(:,3), 5, colorMap_X1)
scatter3(v2_shift(:,1), v2_shift(:,2), v2_shift(:,3), 5, colorMap_X2)

for i = 1:size(X1_i,2)
    ind1 = X1_i(i);
    ind2 = X2_i(i);
    line([v1(ind1,1) v2_shift(ind2,1)],[v1(ind1,2) v2_shift(ind2,2)],[v1(ind1,3) v2_shift(ind2,3)], 'Color', rgb_mat(i,:));
end

hold off
%%

writeVTK('Mesh1.vtk',v1,t1,colorMap_X1);
writeVTK('Mesh2.vtk',v2,t2,colorMap_X2);
writeVTK('Mesh3.vtk',v3,t3,colorMap_X3);
writeVTK('Mesh4.vtk',v4,t4,colorMap_X4);