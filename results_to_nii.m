datadir = 'D:\MRI\SoundPicture\data\raw';
MASK_ORIG_O = struct('subject',num2cell(1:23)', 'filename', {
    fullfile(datadir,'s02_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'03_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'04_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'05_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'06_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'07_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'08_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'09_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'10_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'11_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'12_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'13_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'14_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'15_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'16_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'17_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'18_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'19_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'20_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'21_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'22_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'23_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'24_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
}, 'nii', []);
addpath('dependencies\nifti\');
for i = 1:numel(MASK_ORIG_O)
    MASK_ORIG_O(i).nii = load_nii(MASK_ORIG_O(i).nii);
end