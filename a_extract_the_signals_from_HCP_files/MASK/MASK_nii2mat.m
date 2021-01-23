clc; clear; close all
%HY_96
mkdir('HY_96')
[HY_96_MASK] = y_ReadAll('HarvardOxford-cort-maxprob-thr25-2mm_YCG.nii');% this function from the DPABI used to read the NIFTI file
ROI_label = unique(HY_96_MASK);
for ROI = 1:96
    HY_96 = zeros(91, 109, 91);
    I = 0;
    for x = 1:91
        for y = 1:109
            for z = 1:91
                if HY_96_MASK(x,y,z) == ROI_label(ROI+1)
                    HY_96(x,y,z) = 1;
                    I = I + 1;
                end
            end
        end
    end
    information_HY96MASK(ROI,1).original_label = ROI_label(ROI + 1);
    information_HY96MASK(ROI,1).label = ROI;
    information_HY96MASK(ROI,1).voxcel_number = I;
    save(['HY_96\HY_96_ROI_', num2str(ROI), '.mat'], 'HY_96');
    save('information_HY96MASK.mat', 'information_HY96MASK');
end
%BN_246
mkdir('BN_246')
[BN_246_MASK] = y_ReadAll('BN_Atlas_246_2mm.nii');% this function from the DPABI used to read the NIFTI file
ROI_label = unique(BN_246_MASK);
for ROI = 1:246
    BN_246 = zeros(91, 109, 91);
    I = 0;
    for x = 1:91
        for y = 1:109
            for z = 1:91
                if BN_246_MASK(x, y, z) == ROI_label(ROI + 1)
                    BN_246(x, y, z) = 1;
                    I = I + 1;
                end
            end
        end
    end
    information_BN246MASK(ROI,1).original_label = ROI_label(ROI + 1);
    information_BN246MASK(ROI,1).label = ROI;
    information_BN246MASK(ROI,1).voxcel_number = I;
    save('information_BN246MASK.mat', 'information_BN246MASK');
    save(['BN_246\BN_246_ROI_', num2str(ROI), '.mat'], 'BN_246');
end
%Z-1024
mkdir('Z_1024')
[Z_1024_MASK] = y_ReadAll('Zalesky_1024_parcellated_uniform.nii');% this function from the DPABI used to read the NIFTI file
ROI_label = unique(Z_1024_MASK);
for ROI = 1:1024
    Z_1024 = zeros(91,109,91);
    I = 0;
    for x = 1:91
        for y = 1:109
            for z = 1:91
                if Z_1024_MASK(x,y,z) == ROI_label(ROI + 1)
                    Z_1024(x,y,z) = 1;
                    I = I+1;
                end
            end
        end
    end
    information_Z1024MASK(ROI,1).original_label = ROI_label(ROI + 1);
    information_Z1024MASK(ROI,1).label = ROI;
    information_Z1024MASK(ROI,1).voxcel_number = I;
    save('information_Z1024MASK.mat', 'information_Z1024MASK');
    save(['Z_1024\Z_1024_ROI_', num2str(ROI), '.mat'], 'Z_1024');
end