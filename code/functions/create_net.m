% coord: matrix with 3 colmuns with each row containing MNI coordinates, or a coordinates text file with 4 columns (last as labels)
% NetFlags: meta-structure
% sub: which subject to use as reference (defaults to 1)
% V: reference spm_vol, ideally from a subject
% nifti: string containing nifti file to save with the spheroid ROIs mask, if empty then nifti not saved
% sphere_radius: a number indicating radius of the sphere (if not provided then looks in NetFlags)
% sphere_boundary: a string with how to define the boundaty, options: '<', '<=' (if not provided then looks in NetFlags)
% global_mask_to_use: a string with mask to use (if not provided then looks in NetFlags, if not found then uses 'None')
function create_net(node_size,use_mask,coord,V,nifti)

sphere_radius = node_size;
sphere_boundary = '<=';


if ischar(coord)
    nets = strsplit(coord,'+');
    allfiles = false(1,length(nets));
    for i=1:length(nets)
        if exist(nets{i}, 'file')==2
            allfiles(i) = true;
        end
    end
    if ~all(allfiles)
        nets = {coord};
    end
    
    coord  = [];
    for i=1:length(nets)
        [X Y Z label] = textread(nets{i}, '%f %f %f %s');
        coord = [coord; [X Y Z]];
        %fprintf('coord: %s %d\n', nets{i}, length(X))
    end
end

nROI = size(coord,1);


mask = [];

switch use_mask
    case 'CAT12_02'
        [mask, maskfile, VM] = get_mask(use_mask, '/home/mgell/Work/To_juseless/fingerprint/MaskenEtc/');
        %fprintf('using global mask: %s\n', maskfile)
    case 'FSL025'
        [mask, maskfile, VM] = get_mask(use_mask, '/home/mgell/Work/To_juseless/fingerprint/MaskenEtc/');
        %fprintf('using global mask: %s\n', maskfile)
end

if isempty(V)
    try
        %V = read_data_sub_4D(NetFlags.preproc_data, NetFlags.sub{sub}, NetFlags.SubDir{sub}, NetFlags.datadir, NetFlags.data{1});
        file4d = sub_files(NetFlags.preproc_data, NetFlags.sub{sub}, NetFlags.SubDir{sub}, NetFlags.data{1});
        V = spm_vol(fullfile(NetFlags.datadir, file4d));
    catch ME
        disp('data could not be read')
        throw(ME)
    end
end

% get all indices
clearvars XYZ;
tmp = ones(V(1).dim(1:3)); ind = find(tmp);
[XYZ(:,1), XYZ(:,2), XYZ(:,3)] = ind2sub(size(tmp),ind);
XYZ = XYZ';
XYZmm = V(1).mat * [XYZ; ones(1,size(XYZ,2))];
XYZmm = XYZmm(1:3,:);

%fprintf('using spheres %s %dmm\n', sphere_boundary, sphere_radius)
% get the spheres, this will include all gm, wm and csf
Ind = []; VOI = [];
for ana = 1:nROI
    nodeXYZmm = coord(ana,:); % this must be already in mm
    D = pdist2(nodeXYZmm, XYZmm'); % dist of each voxel to the node center
    assert(length(D)==size(XYZmm,2))
    % get indices of voxels inside the sphere
    switch sphere_boundary
        case '<='
            sphere_vox = find(D<=sphere_radius);
        case '<'
            sphere_vox = find(D<sphere_radius);
        otherwise
            error(['Unknwon sphere boundary: ' sphere_boundary])
    end
    Ind = [Ind sphere_vox];
    VOI = [VOI repmat(ana,1,length(sphere_vox))];
end
assert(length(Ind)==length(VOI))

if ~isempty(mask)
    [maskXYZ(:,1), maskXYZ(:,2), maskXYZ(:,3)] = ind2sub(VM.dim, find(mask > 0));
    maskXYZ = maskXYZ';
    maskXYZmm = VM.mat * [maskXYZ; ones(1, size(maskXYZ,2))];
    maskXYZmm = maskXYZmm(1:3,:);
    % match sphere voxels with the mask if used
    % get the min distance to mask in mm and subset the Ind and VOI with that
    D_min_mask = pdist2(maskXYZmm', XYZmm(:,Ind)','euclidean','Smallest',1);
    assert(length(D_min_mask)==length(Ind)) % we must find all coordinates?
    matched_vox = D_min_mask==0;
    Ind = Ind(matched_vox);
    VOI = VOI(matched_vox);
end

assert(length(Ind)==length(VOI))
assert(size(XYZ,2)>=max(Ind))

if ~isempty(nifti)
    assert(ischar(nifti))
    V = V(1);
    Vdat = zeros(V.dim(1:3));
    uVOI = unique(VOI);
    for i=1:length(uVOI)
        Vdat(Ind(VOI==uVOI(i))) = uVOI(i);
    end
    V.fname = nifti;
    V.private.dat.fname = nifti;
    [~] = spm_write_vol(V, Vdat);
end

end



