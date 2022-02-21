
%%$ CREATE RANDOM NETOWORKS USNG PERMUTATION
% Requirements:
%   1. spm12 and NIFTI toolboxes
%   2. Template network: netowrk used as template for ranom networks
%   3. Atlas file of nodes if not using whole brain
%   3. One subject file with same resolution and orientation
%       as random network nifti to view results, placed in 'sub_files' dir

clear

% paths to spm and Nifti toolbox
addpath('/home/mgell/Matlab_toolboxes/spm12');
addpath('/home/mgell/Matlab_toolboxes/NIFTI_toolbox-master');

% directory with project
wd = '/home/mgell/Work/random_nets';
addpath(genpath(wd));
cd(wd);


%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_net = 'networks/MotorHeckner_VOIs.txt'
Atlas = 'networks/power_nodes.node';

whole_brain = true % True for hole brain, False for sampling from nodes.

node_size = 5 % sphere size
use_mask = 'FSL025' % '' or 'CAT12_02' or 'FSL025'
nperm = 10 % n of random nets
community = ''; % whole network or specific community? Need to provide community assignment. 
                % Use same value as in community assignment. e.g. 1 if only
                % want community called 1. Requires column 5 with community info in 'nodes'.

all_masks = false % True for print all networks as nii. If True set below:
ref_nii = 'sub_files/rfMRI_REST1_AP_hp0_clean.nii.gz' % used as a reference for creating nii with nodes
                % dimension and orientation information is taken from this file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load stuff
% Template network
temp_net = readtable(fullfile(wd, temp_net));

% Atlas to sample from - needs to be xyz coordiantes + name/id in last column
if whole_brain == false
    nodes = readtable(fullfile(wd, Atlas), 'FileType','text');
    node_id = nodes(:,4);
end

% Output dir
outputdir = fullfile(wd, 'res');

% Load reference nifti
V = spm_vol(fullfile(wd, ref_nii));

% subsample network to community if requested
if ~isempty(community)
    communities = readtable('/home/mgell/Work/FC/hcp/text_files/res/replic_5mm_net_communities.csv','ReadVariableNames',false);
    temp_net(:,5) = communities;
    rows = temp_net.Var5 == community;
    temp_net = temp_net(rows,:);
end

node_n = size(temp_net,1);
all_y = table(zeros(node_n,1));



%%% Do the thing
for perm_i = 1:nperm

    fprintf('permutation: %i\n', perm_i) % permutation number
    
    duplicate = 1;
    while duplicate == 1 
        % create random net
        if whole_brain == false
            [y, nodes_i] = nodes_just_right(nodes,node_id,node_n,temp_net);
        else
            [y, nodes_i] = nodes_just_right_whole_brain(node_n,temp_net);
            nodes_i = array2table([nodes_i, y],'VariableNames',{'x','y','z','id'});
        end
        
        % check if random network already exists
        duplicate = check_duplicates(y,all_y); 
    end
    
    all_y(:,perm_i) = array2table(y); % save all ys so we can avoid duplicates of random samples
    
    % save random network
    coords_i = fullfile(outputdir, ['coords' num2str(perm_i) '.txt']);
    if all_masks == true
        net_i = fullfile(outputdir, ['net' num2str(perm_i) '.nii']);
    else
        net_i = fullfile(outputdir, 'net.nii');
    end

    % Save
    writetable(nodes_i,coords_i,'WriteVariableNames',false,'Delimiter','tab')
    create_net(node_size,use_mask,coords_i,V,net_i);
end

