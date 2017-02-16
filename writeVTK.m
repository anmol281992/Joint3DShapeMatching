function writeVTK(filename, vertices, faces, func_vec)

currdir = pwd;
filename = strcat(currdir,'\VTKFiles\',filename);
if size(vertices, 1) < size(vertices, 2)
   vertices = vertices';
end

if size(faces, 1) < size(faces, 2)
    faces = faces';
end

numpts = size(vertices, 1);
numfaces = size(faces, 1);

nl = sprintf('\n');

ofid = fopen(filename, 'wt');

fwrite(ofid, ['# vtk DataFile Version3.0' nl 'vtk output' nl 'ASCII' nl ...
    'DATASET POLYDATA' nl 'POINTS ' num2str(numpts) ' float' nl]);

for i = 1:numpts
    fwrite(ofid, [num2str(vertices(i,1)) ' ' num2str(vertices(i,2)) ' ' ...
        num2str(vertices(i,3)) nl]);
end

fwrite(ofid, ['POLYGONS ' num2str(numfaces) ' ' num2str(numfaces*4) nl]);

faces = faces-1;
for i = 1:numfaces
    fwrite(ofid, ['3 ' num2str(faces(i,1)) ' ' num2str(faces(i,2)) ' ' ...
        num2str(faces(i,3)) nl]);
end

fwrite(ofid, nl);
fwrite(ofid, ['POINT_DATA ' num2str(numpts) nl ...
    'SCALARS distance_from float' nl 'LOOKUP_TABLE default' nl]);
for i = 1:numpts
    fwrite(ofid, [num2str(func_vec(i)) nl]);
end

fclose(ofid);