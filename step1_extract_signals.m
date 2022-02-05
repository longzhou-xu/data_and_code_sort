%% step1: extract the ROI signals from the HCP data
clc;
clear;
close all;

%% 1. set the file path
% set the HCP rfMRI data storage path

HCP_path = 'the folder containing the folder';
subj_path = dir(fullfile(HCP_path,'*'));
rsfMRI_fname = 'rfMRI_REST1_LR_hp2000_clean.nii.gz';

%% 2. use the atlas to extract the ROI signals
for SUB = 1:length(subj_path)
    
    % show the subject in command line
    display(['subject ID is ',subj_path(SUB, 1).name])
    ROIsignals(SUB).subj_ID = str2double(subj_path(SUB, 1).name);
    
    % import the rfMRI data into MATLAB workspace
    rfMRI_path = fullfile(HCP_path, subj_path(SUB).name, rsfMRI_fname);
    % load the data, "y_ReadAll" from the DPABI package;
    rfMRI = y_ReadAll(rfMRI_path);
    
    % using the Harvard-oxford cortical brain atlas to extract the 96 ROI signals
    for ROI = 1:96
        
        % load the Harvard-oxford atlas
        ROIname = ['HY_96_ROI_', num2str(ROI), '.mat'];
        load(fullfile('atlas','HY_96',ROIname));
        
        % extract the voxel signals
        clear voxels_ROI
        ROIdata = rfMRI(:, :, :, 1) .* HY_96;
        ROIindex = f_find_3d(ROIdata);
        voxels_ROI(ROI).coor = ROIindex;
        for I = 1:length(ROIindex)
            voxels_ROI(ROI).signals(I, :) = zscore(rfMRI(ROIindex(I, 1), ROIindex(I, 2), ROIindex(I, 3), :));
        end
        
        % save the voxel signals
        mkdir(fullfile('step_1','voxel_signals','HY_96', subj_path(SUB,1).name)); % folder for saving the voxel signals
        save(fullfile('step_1','voxel_signals','HY_96', subj_path(SUB,1).name,'voxel_signals.mat'),'voxels_ROI','-v7.3');
        
        % extract the ROI signals
        ROIsignals_HY(SUB).ROIsignals(:, ROI) = mean(voxels_ROI(ROI).signals, 1);
        
    end
   %save the ROI signals
    mkdir(fullfile('step_1','ROI_signals'));
    save(fullfile('step_1','ROI_signals','ROIsignals.mat'), 'results','-v7.3');
    
    % using the BN-246 brain atlas to extract the BN-246 ROI signals
    for ROI = 1:246
        
        % load the Harvard-oxford atlas
        ROIname = ['BN_246_ROI_', num2str(ROI), '.mat'];
        load(fullfile('atlas','BN_246',ROIname));
        
        % extract the voxel signals
        clear voxels_ROI
        ROIdata = rfMRI(:, :, :, 1) .* BN_246;
        ROIindex = f_find_3d(ROIdata);
        voxels_ROI(ROI).coor = ROIindex;
        for I = 1:length(ROIindex)
            voxels_ROI(ROI).signals(I, :) = zscore(rfMRI(ROIindex(I, 1), ROIindex(I, 2), ROIindex(I, 3), :));
        end
        
        % save the voxel signals
        mkdir(fullfile('step_1','voxel_signals','BN_246', subj_path(SUB,1).name)); % folder for saving the voxel signals
        save(fullfile('step_1','voxel_signals','BN_246', subj_path(SUB,1).name,'voxel_signals.mat'),'voxels_ROI','-v7.3');
        
        % extract the ROI signals
        ROIsignals(SUB).ROIsignals_BN(:, ROI) = mean(voxels_ROI(ROI).signals, 1);
        
    end
   %save the ROI signals
    mkdir(fullfile('step_1','ROI_signals'));
    save(fullfile('step_1','ROI_signals','ROIsignals.mat'), 'results','-v7.3');
    
    % using the Z-1024 brain atlas to extract the Z-1024 ROI signals
    for ROI = 1:1024
        
        % load the Harvard-oxford atlas
        ROIname = ['Z_1024_ROI_', num2str(ROI), '.mat'];
        load(fullfile('atlas','Z_1024',ROIname));
        
        % extract the voxel signals
        clear voxels_ROI
        ROIdata = rfMRI(:, :, :, 1) .* Z_1024;
        ROIindex = f_find_3d(ROIdata);
        voxels_ROI(ROI).coor = ROIindex;
        for I = 1:length(ROIindex)
            voxels_ROI(ROI).signals(I, :) = zscore(rfMRI(ROIindex(I, 1), ROIindex(I, 2), ROIindex(I, 3), :));
        end
        
        % save the voxel signals
        mkdir(fullfile('step_1','voxel_signals','Z_1024', subj_path(SUB,1).name)); % folder for saving the voxel signals
        save(fullfile('step_1','voxel_signals','Z_1024', subj_path(SUB,1).name,'voxel_signals.mat'),'voxels_ROI','-v7.3');
        
        % extract the ROI signals
        ROIsignals(SUB).ROIsignals_Z(:, ROI) = mean(voxels_ROI(ROI).signals, 1);
        
    end
   %save the ROI signals
    mkdir(fullfile('step_1','ROI_signals'));
    save(fullfile('step_1','ROI_signals','ROIsignals.mat'), 'results','-v7.3');
    
end
% save the subject ID
save('step_1\subject.mat','subject');
%% function
function [index] = f_find_3d(data)
%[index] = find_nonzero_3d(data)
% This function is used to find the index of non-zero elements in a 3D matrix.
    Size=size(data);
    index=[];
    I=0;
    for x=1:Size(1)
        for y=1:Size(2)
            for z=1:Size(3)
                if data(x,y,z)~=0
                    I=I+1;
                    index(I,1)=x;
                    index(I,2)=y;
                    index(I,3)=z;
                end
            end
        end
    end
end