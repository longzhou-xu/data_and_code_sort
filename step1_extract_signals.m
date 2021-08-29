%% step1: extract the ROI signal from the HCP data

%% setting
clc;clear;close all;
% environment: DPABI,spm

%% 1. set the file path
% set the HCP rfMRI data storage path
% the data path: subject_folder_path\subject_folder\rfMRI_filename
subject_folder_path = 'HCP_RS_LR';
subject_folder = dir(subject_folder_path);
subject_folder = subject_folder(3:end);
rfMRI_filename = 'rfMRI_REST1_LR_hp2000_clean.nii.gz';

%% 2. use the atlas to extract the ROI signals
for SUB = 1:length(subject_folder)
    
    % show the subject in command line
    display(['subject ID is ',subject_folder(SUB, 1).name])
    subject(SUB,1) = str2double(subject_folder(SUB, 1).name);
    
    % load the rfMRI data
    % the data path
    clear rfMRIdatapath
    rfMRIdatapath = [subject_folder_path,'\',subject_folder(SUB).name,'\',rfMRI_filename];
    % laod the data, note that the y_ReadAll from the DPABI package;
    rfMRIdata = y_ReadAll(rfMRIdatapath);
    
    % using the Harvard-oxford cortical brain atlas to extract the 96 ROI signals
    for ROI = 1:96
        
        %extract the voxcel signals
        % setting
        clear RS_HY_96_voxceL_signals index data_aftermask
        % load the Harvard-oxford atlas
        load(['step_1\MASK\HY_96\HY_96_ROI_', num2str(ROI), '.mat']);
        % extrac the voxcel signal
        data_aftermask = rfMRIdata(:, :, :, 1) .* HY_96;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_HY_96_voxceL_signals(I, :) = rfMRIdata(index(I, 1), index(I, 2), index(I, 3), :);
        end
        % save the voxcel signals
        mkdir(['step_1\voxcel_siganls\HY_96\', subject_folder(SUB,1).name,...
            '\ROI',num2str(ROI)]);
        save(['step_1\voxcel_siganls\HY_96\', subject_folder(SUB,1).name,...
            '\ROI',num2str(ROI),'\voxcel_signals.mat'],'RS_HY_96_voxceL_signals');
        save(['step_1\voxcel_siganls\HY_96\', subject_folder(SUB,1).name,...
            '\ROI',num2str(ROI),'\index_HY96.mat'],'index');
        
        % extract the ROI signals
        rest_HY96_ROI(:, ROI) = mean(RS_HY_96_voxceL_signals, 1);
    end
   %save the ROI signals
    mkdir('step_1\ROI_signals\HY_96');
    save(['step_1\ROI_signals\HY_96\sub', num2str(SUB), '.mat'], 'rest_HY96_ROI');
    
    %extract the signal using the BN-246 atlas
    for ROI = 1:246
        
        % extract the voxcel signals
        % setting
        clear RS_BN_246_voxceL_signals index data_aftermask
        %laod the mask
        load(['step_1\MASK\BN_246\BN_246_ROI_', num2str(ROI), '.mat']);
        %extract the voxcel siganls
        data_aftermask = rfMRIdata(:, :, :, 1) .* BN_246;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_BN_246_voxceL_signals(I, :) = rfMRIdata(index(I, 1), index(I, 2), index(I, 3), :);
        end
        %save the voxcel signals
        mkdir(['step_1\voxcel_siganls\BN_246\', subject_folder(SUB,1).name, '\ROI',num2str(ROI)]);
        save(['step_1\voxcel_siganls\BN_246\', subject_folder(SUB,1).name, '\ROI',num2str(ROI), '\voxcel_signals.mat'], 'RS_BN_246_voxceL_signals');
        save(['step_1\voxcel_siganls\BN_246\', subject_folder(SUB,1).name, '\ROI',num2str(ROI), '\index_BN246.mat'], 'index');
        
        %extract the ROI signals
        ROI_BN246_RS(:, ROI) = mean(RS_BN_246_voxceL_signals, 1);
    end
    % save the ROI signals
    mkdir('step_1\ROI_signals\BN_246');
    save(['step_1\ROI_signals\BN_246\sub', num2str(SUB), '.mat'], 'ROI_BN246_RS');
    
    %extract the signa using the Z_1024 atlas
    for ROI = 1:1024
        
        % extract the voxcel siganls
        % setting
        clear RS_Z_1024_voxceL_signals index data_aftermask
        % load the mask
        load(['step_1\MASK\Z_1024\Z_1024_ROI_', num2str(ROI), '.mat']);
        % extract the siganls
        data_aftermask = rfMRIdata(:, :, :, 1) .* Z_1024;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_Z_1024_voxceL_signals(I, :) = rfMRIdata(index(I, 1), index(I, 2), index(I, 3), :);
        end
        %save the voxcel signals
        mkdir(['step_1\voxcel_siganls\Z_1024\', subject_folder(SUB,1).name, '\ROI',num2str(ROI)]);
        save(['step_1\voxcel_siganls\Z_1024\', subject_folder(SUB,1).name, '\ROI',num2str(ROI), '\voxcel_signals.mat'], 'RS_Z_1024_voxceL_signals');
        save(['step_1\voxcel_siganls\Z_1024\', subject_folder(SUB,1).name, '\ROI',num2str(ROI), '\index_Z1024.mat'], 'index');
        
        %extract the ROI signals
        ROI_Z1024_RS(:, ROI) = mean(RS_Z_1024_voxceL_signals, 1);
    end
    % save the ROi signal
    mkdir('step_1\ROI_signals\Z_1024');
    save(['step_1\ROI_signals\Z_1024\sub', num2str(SUB), '.mat'], 'rest_Z1024_ROI');
end
% save the subject ID
save('step_1\subject.mat','subject');
%% function
function [index] = find_nonzero_3d(data)
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