rootDir = pwd();
addpath(fullfile(rootDir, 'dependencies', 'nifti'));
datadir = '/Users/Chris/MRI/Manchester/data/raw';
MASK_ORIG_O = struct('subject',num2cell(1:23)', 'filename', {
    fullfile(datadir,'s02_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s03_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s04_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s05_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s06_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s07_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s08_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s09_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s10_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s11_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s12_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s13_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s14_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s15_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s16_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s17_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s18_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s19_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s20_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s21_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s22_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s23_leftyes','mask','nS_c1_mask_nocerebellum_O.nii')
    fullfile(datadir,'s24_rightyes','mask','nS_c1_mask_nocerebellum_O.nii')
}, 'hdr', []);
for i = 1:numel(MASK_ORIG_O)
    MASK_ORIG_O(i).hdr = load_nii_hdr(MASK_ORIG_O(i).filename);
end