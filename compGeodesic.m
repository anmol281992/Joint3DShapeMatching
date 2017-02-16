function compGeodesic(path)

filelist = dir(path);
L = 1;
for i = 4:numel(filelist)
    b = filelist(i).name
    if(strcmp(b(end-3:end),'.off'))
        %read mesh from the off file
        fid=fopen(b);
        fgetl(fid);
        nos = fscanf(fid, '%d %d  %d', [3 1]);
        nopts = nos(1);
        notrg = nos(2);
        coord = fscanf(fid, '%g %g  %g', [3 nopts]);
        coord = coord';
        triang=fscanf(fid, '%d %d %d %d',[4 notrg]);
        triang=triang';
        triang=triang(:,2:4)+1; %%we have added 1 because the vertex indices start from 0 in vtk format
        fclose(fid);
        
        % compute geodesic
        numpts = size(coord,1);
        Edist = pdist2(coord, coord, 'euclidean');
    
        G = zeros(numpts,numpts);
        for j = 1:size(triang, 1)
            p1 = triang(j,1);
            p2 = triang(j,2);
            p3 = triang(j,3);
    
            G(p1,p2)= Edist(p1,p2);
            G(p2,p1)= Edist(p2,p1);
            G(p1,p3)= Edist(p1,p3);
            G(p3,p1)= Edist(p3,p1);
            G(p2,p3)= Edist(p2,p3);
            G(p3,p2)= Edist(p3,p2);
        end
        geodesic(L).Gdist = graphallshortestpaths(sparse(G));
        L = L+1;
    end
end

save('geodesic','geodesic');