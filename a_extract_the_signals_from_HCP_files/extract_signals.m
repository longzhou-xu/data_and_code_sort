%this code translate the HCP resting data to the mat data
clc; clear; close all
file_folder = dir('HCP_RS_LR'); %the path of the RS-fMRI data from HCP
file_folder = file_folder(3:end);
for S = 1:length(file_folder)
    file_folder(S, 1).name
    data = y_ReadAll(['HCP_RS_LR\', ...
        file_folder(S).name, '\rfMRI_REST1_LR_hp2000_clean.nii.gz']); % this function from the DPABI used to read the NIFTI file
    % HY_96
    for roi = 1:96
        clear RS_HY_96_voxceL_signals index data_aftermask
        load(['MASK\HY_96\HY_96_ROI_', num2str(roi), '.mat']);
        mkdir(['voxcel_siganls\HY_96\', file_folder(S,1).name, '\ROI',num2str(roi)]);
        data_aftermask = data(:, :, :, 1) .* HY_96;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_HY_96_voxceL_signals(I, :) = data(index(I, 1), index(I, 2), index(I, 3), :);
        end
        save(['voxcel_siganls\HY_96\', file_folder(S,1).name, '\ROI',num2str(roi),'\voxcel_signals.mat'],'RS_HY_96_voxceL_signals');
        save(['voxcel_siganls\HY_96\', file_folder(S,1).name, '\ROI',num2str(roi),'\index_HY96.mat'],'index');
        rest_HY96_ROI(:, roi) = mean(RS_HY_96_voxceL_signals, 1);
    end
    mkdir('ROI_signals\HY_96');
    save(['ROI_signals\HY_96\sub', num2str(S), '.mat'], 'rest_HY96_ROI');
    
    %BN_246
    for roi = 1:246
        clear RS_BN_246_voxceL_signals index data_aftermask
        load(['MASK\BN_246\BN_246_ROI_', num2str(roi), '.mat']);
        mkdir(['voxcel_siganls\BN_246\', file_folder(S,1).name, '\ROI',num2str(roi)]);
        data_aftermask = data(:, :, :, 1) .* BN_246;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_BN_246_voxceL_signals(I, :) = data(index(I, 1), index(I, 2), index(I, 3), :);
        end
        save(['voxcel_siganls\BN_246\', file_folder(S,1).name, '\ROI',num2str(roi), '\voxcel_signals.mat'], 'RS_BN_246_voxceL_signals');
        save(['voxcel_siganls\BN_246\', file_folder(S,1).name, '\ROI',num2str(roi), '\index_BN246.mat'], 'index');
        rest_BN246_ROI(:, roi) = mean(RS_BN_246_voxceL_signals, 1);
    end
    mkdir('ROI_signals\BN_246');
    save(['ROI_signals\BN_246\sub', num2str(S), '.mat'], 'rest_BN246_ROI');
    
    %Z_1024
    for roi = 1:1024
        clear RS_Z_1024_voxceL_signals index data_aftermask
        load(['MASK\Z_1024\Z_1024_ROI_', num2str(roi), '.mat']);
        mkdir(['voxcel_siganls\Z_1024\', file_folder(S,1).name, '\ROI',num2str(roi)]);
        data_aftermask = data(:, :, :, 1) .* Z_1024;
        index = find_nonzero_3d(data_aftermask);
        for I = 1:length(index)
            RS_Z_1024_voxceL_signals(I, :) = data(index(I, 1), index(I, 2), index(I, 3), :);
        end
        save(['voxcel_siganls\Z_1024\', file_folder(S,1).name, '\ROI',num2str(roi), '\voxcel_signals.mat'], 'RS_Z_1024_voxceL_signals');
        save(['voxcel_siganls\Z_1024\', file_folder(S,1).name, '\ROI',num2str(roi), '\index_Z1024.mat'], 'index');
        rest_Z1024_ROI(:, roi) = mean(RS_Z_1024_voxceL_signals, 1);
    end
    mkdir('ROI_signals\Z_1024');
    save(['ROI_signals\Z_1024\sub', num2str(S), '.mat'], 'rest_Z1024_ROI');
end