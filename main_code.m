% place all the mesh files in the Data directory

if(exist('./Data/hkslist.mat', 'file') ~= 2)
    compHKS('.\Data');
end

if(exist('./Data/geodesic.mat', 'file') ~=2)
    compGeodesic('./Data');
end

datadir = strcat(pwd, '/Data');
filelist = dir(datadir);

L = 1;
for i = 3:numel(filelist)
    b = filelist(i).name;
    filepath = strcat(datadir,'/',b);
    if(strcmp(filepath(end-3:end),'.off'))
        fid=fopen(filepath);
        fgetl(fid);
        nos = fscanf(fid, '%d %d  %d', [3 1]);
        nopts = nos(1);
        notrg = nos(2);
        coord = fscanf(fid, '%g %g  %g', [3 nopts]);
        p(L).name = b(1:end-4);
        p(L).coord = coord';
        triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
        triang=triang';
        p(L).triang=triang(:,2:4)+1;
    
        L = L+1;
        fclose(fid);
    end
end

L = L-1;
points = [1:125:12500];

load('.\Data\hkslist');
load('.\Data\geodesic');

%% Create S matrices

s = zeros(100,101,L);
S = zeros(100,100,L);
hks_i = zeros(12500,101);

for i = 1:L
    hks_i = P(i).hks;
    s(:,:,i) = hks_i(points,:);
end

Srel2 = [];

for n = 1:L
    Srel1 = [];
    for m = 1:L
        s1 = s(:,:,m);
        Srel = [];
        for i = 1:100
            S  = sum((s1 - repmat(s1(i,:),100,1)).^2,2);
            Srel = [Srel; S'];
        end
        Srel1 = [Srel1 Srel];
    end
    Srel2 = [Srel2; Srel1];
end

%% Run Match ALS
alpha = 0.2;
dimGroup = ones(1,L) * 100;
Wrel = -Srel2 + alpha*(ones(L*100,L*100));
W = sparse(Wrel);
[X,info,A] = mmatch_CVX_ALS(W,dimGroup);
X = full(X);

%% Distrotion Term
addpath gm-toolbox

Xdistort = zeros(L*100, L*100);
for i = 1:L
    for j = 1:L
        iStrt = (i-1)*100 + 1;
        jStrt = (j-1)*100 + 1;
        Xij = X(iStrt:iStrt+99, jStrt:jStrt+99);
        Sij = Srel2(iStrt:iStrt+99, jStrt:jStrt+99);
        [N,M] = find(Xij == 1);
        match = [N';M'];
        n = length(N);
        sim = zeros(1,n);
        for k = 1:n
            sim(k) = Sij(N(k),M(k));
        end
        
        %% solve graph matching
        wEdge = 10;
        if wEdge <= 0
            Xraw = double(sim);
        else
            % construct graphs
            fprintf('- constructing graphs ...\n');
            [uniq_feat1,~,new_feat1] = unique(match(1,:));
            [uniq_feat2,~,new_feat2] = unique(match(2,:));
            cand_matchlist_uniq = [new_feat1';new_feat2'];
            xyi = (p(i).coord)';
            xyj = (p(j).coord)';
            edgeAttr1 = computeEdgeAttr(xyi(:,uniq_feat1));
            edgeAttr2 = computeEdgeAttr(xyj(:,uniq_feat2));
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
        
        xijdistort = greedyMatch(match,Xraw);
        Xijfind = find(xijdistort);
        match1 = match(:,Xijfind);

        Xijdistort = zeros(100,100);
        [n1,m1] = size(match1);
        for k = 1:m1
            Xijdistort(match1(1,k),match1(2,k)) = 1;
        end

        Xdistort(iStrt:iStrt+99, jStrt:jStrt+99) = Xijdistort;
    end
end


%% perform color matching for possible pairs

for i=1:L
    for j=i+1:L
        iStrt = (i-1)*100 + 1;
        jStrt = (j-1)*100 + 1;
        Xdistortij = Xdistort(iStrt:iStrt+99, jStrt:jStrt+99);
        v1 = p(i).coord;    t1 = p(i).triang;
        v2 = p(j).coord;    t2 = p(j).triang;
        Gdist1 = geodesic(i).Gdist;
        Gdist2 = geodesic(j).Gdist;
        n = size(v1,1);
        [Xi, Xj] = find(Xdistortij);
        Xi = points(Xi);
        Xj = points(Xj);
        [colormapi, colormapj] = assignColors(n, Gdist1, Gdist2, Xi, Xj);
        
        % write VTK files for models and scatter plot.
        filenamei = sprintf('%s_X%d%d.vtk', p(i).name, i, j);
        filenamej = sprintf('%s_X%d%d.vtk', p(j).name, i, j);
        writeVTK(filenamei, v1, t1, colormapi);
        writeVTK(filenamej, v2, t2, colormapj);
        
        figure()
        subplot(1,2,1)
        scatter3(v1(:,1), v1(:,2), v1(:,3), 3, colormapi)
        title(p(i).name);

        subplot(1,2,2)
        scatter3(v2(:,1), v2(:,2), v2(:,3), 3, colormapj)
        title(p(j).name);
        
        %% create file for evaluation
        
        [k, l] = find(Xdistortij);
        list = [points(k); points(l)];
        currdir = pwd;
        filename = sprintf('%s_%s_X%d%d.cor', p(i).name, p(j).name, i, j);
        filename = strcat(currdir, '\Evaluation\CORFiles\',filename)
        ofid = fopen(filename, 'wt');
        for count = 1:length(list)
            fprintf(ofid,'%d\n',list(count));
        end
    end
end