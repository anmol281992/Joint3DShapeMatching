clear all;
close all;
filelist = dir('/Users/compume/Documents/MATLAB/computer_vision_project/release');

% all noise in the mesh was added in this code error induced preventing it
% from running remove error k = 
L = 1;
for i = 4:numel(filelist)
a = filelist(i).name;
b = '/Users/compume/Documents/MATLAB/computer_vision_project/release/';
b = strcat(b,a);
%p =struct([]);
if(strcmp(b(end-3:end),'.off'))
    fid=fopen(b);
    fgetl(fid);
    nos = fscanf(fid, '%d %d  %d', [3 1]);
    p(L).nopts = nos(1);
    p(L).notrg = nos(2);
    coord = fscanf(fid, '%g %g  %g', [3 p(L).nopts]);
    p(L).field = a;
    p(L).coord = coord';
    triang=fscanf(fid, '%d %d %d %d',[4 p(L).notrg]);
    triang=triang';
    p(L).triang=triang(:,2:4)+1;
    
    L = L+1;
    fclose(fid);
end
%%we have added 1 because the vertex indices start from 0 in vtk format


end
%load('centaur_HKS.mat');
figure(1);
hold on
plot3(p(1).coord(:,1),p(1).coord(:,2),p(1).coord(:,3),'g.');

plot3(p(1).coord(1:200:3400,1),p(1).coord(1:200:3400,2),p(1).coord(1:200:3400,3),'r.');
hold off
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
       S2(:,i,3)  = (sum((s2 - repmat(s3(i,:),100,1)).^2,2));
       S2(:,i,4)  = (sum((s2 - repmat(s4(i,:),100,1)).^2,2));
    end
    for i = 1:100
       S3(:,i,1)  = (sum((s3 - repmat(s1(i,:),100,1)).^2,2));
       S3(:,i,2)  = (sum((s3 - repmat(s2(i,:),100,1)).^2,2));
       S3(:,i,3)  = (sum((s3 - repmat(s3(i,:),100,1)).^2,2));
       S3(:,i,4)  = (sum((s3 - repmat(s4(i,:),100,1)).^2,2));
    end
    for i = 1:100
       S4(:,i,1)  = (sum((s4 - repmat(s1(i,:),100,1)).^2,2));
       S4(:,i,2)  = (sum((s4 - repmat(s2(i,:),100,1)).^2,2));
       S4(:,i,3)  = (sum((s4 - repmat(s3(i,:),100,1)).^2,2));
       S4(:,i,4)  = (sum((s4 - repmat(s4(i,:),100,1)).^2,2));
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

[k,m] = find(X12 == 1);
%f_vec1 = zeros(1,12500);

% fvec1(points(k(j)):(points(k(j))+123)) = colorvec(k(j));
 %fvec2(points(m(j)):(points(m(j))+123)) = colorvec(k(j));
%end%
%% We have found the assignment 
%% Now we have to match the assignment with color
% % Geodesic Distance
coord1 = p(1).coord;
nopts1 = p(1). nopts;
notrg1 = p(1). notrg;
triang1 = p(1).triang;
        
coord2 = p(2).coord;
nopts2 = p(2). nopts;
triang2 = p(2).triang;
notrg2 = p(2). notrg;
    
% COLORMAP CREATION
LIST = [1:125:12500];
correspond = zeros(1,100);
source = 1;
descriptorgeodesic1 = zeros(12500,12500);
    [v,c] = size(triang1);
 for i = 1:v
        descriptorgeodesic1(triang1(i,1),triang1(i,2)) = sqrt(sum((coord1(triang1(i,2),:)-coord1(triang1(i,1),:)).^2));
        descriptorgeodesic1(triang1(i,1),triang1(i,3)) = sqrt(sum((coord1(triang1(i,3),:)-coord1(triang1(i,1),:)).^2));
        descriptorgeodesic1(triang1(i,2),triang1(i,3)) = sqrt(sum((coord1(triang1(i,3),:)-coord1(triang1(i,2),:)).^2));
        descriptorgeodesic1(triang1(i,2),triang1(i,1)) = descriptorgeodesic1(triang1(i,1),triang1(i,2));
        descriptorgeodesic1(triang1(i,3),triang1(i,1)) = descriptorgeodesic1(triang1(i,1),triang1(i,3));
        descriptorgeodesic1(triang1(i,3),triang1(i,2)) = descriptorgeodesic1(triang1(i,2),triang1(i,3));
 end
 descriptorgeodesic2 = zeros(12500,12500);
    [v,c] = size(triang2);
    for i = 1:v
        descriptorgeodesic2(triang2(i,1),triang2(i,2)) = sqrt(sum((coord2(triang2(i,2),:)-coord2(triang2(i,1),:)).^2));
        descriptorgeodesic2(triang2(i,1),triang2(i,3)) = sqrt(sum((coord2(triang2(i,3),:)-coord2(triang2(i,1),:)).^2));
        descriptorgeodesic2(triang2(i,2),triang2(i,3)) = sqrt(sum((coord2(triang2(i,3),:)-coord2(triang2(i,2),:)).^2));
        descriptorgeodesic2(triang2(i,2),triang2(i,1)) = descriptorgeodesic2(triang2(i,1),triang2(i,2));
        descriptorgeodesic2(triang2(i,3),triang2(i,1)) = descriptorgeodesic2(triang2(i,1),triang2(i,3));
        descriptorgeodesic2(triang2(i,3),triang2(i,2)) = descriptorgeodesic2(triang2(i,2),triang2(i,3));
    end
L1 = zeros(100,12500);
L2 = zeros(100,12500);
for source = 1:1:100
L1(source,:) = graphshortestpath(sparse(descriptorgeodesic1),LIST(source));
source1 = find(X12(source,:)== 1);
correspond(source1) = source;
L2(source1,:) = graphshortestpath(sparse(descriptorgeodesic2),LIST(source));
end
closeness = [1:1:100];
[Min,I] = min(L1,[],1);
f_vec1 = closeness(I);
[Min,I] = min(L2,[],1);
f_vec2 = closeness(correspond(I));

ofid = fopen('man1.vtk','w');
    fprintf(ofid, '# vtk DataFile Version 3.0\n');
    fprintf(ofid,'vtk output\n');
    fprintf(ofid,'ASCII\n');
    fprintf(ofid,'DATASET POLYDATA\n');
    fprintf(ofid,'POINTS %d float\n', nopts);
    fprintf(ofid,'%g %g %g\n', coord');
    fprintf(ofid,'POLYGONS %d %d\n', notrg, 4*notrg);
    fprintf(ofid,'3 %d %d %d\n', triang'-1);
    fprintf(ofid,'\n');
    fprintf(ofid,'POINT_DATA %d\n', nopts);
    fprintf(ofid,'SCALARS distance_from float\n');
    fprintf(ofid,'LOOKUP_TABLE default\n');
    fprintf(ofid,'%g\n',  f_vec1');
    fclose(ofid);
    
    ofid = fopen('man2.vtk','w');
    fprintf(ofid, '# vtk DataFile Version 3.0\n');
    fprintf(ofid,'vtk output\n');
    fprintf(ofid,'ASCII\n');
    fprintf(ofid,'DATASET POLYDATA\n');
    fprintf(ofid,'POINTS %d float\n', nopts1);
    fprintf(ofid,'%g %g %g\n', coord1');
    fprintf(ofid,'POLYGONS %d %d\n', notrg1, 4*notrg1);
    fprintf(ofid,'3 %d %d %d\n', triang1'-1);
    fprintf(ofid,'\n');
    fprintf(ofid,'POINT_DATA %d\n', nopts1);
    fprintf(ofid,'SCALARS distance_from float\n');
    fprintf(ofid,'LOOKUP_TABLE default\n');
    fprintf(ofid,'%g\n',  f_vec2');
    fclose(ofid);




