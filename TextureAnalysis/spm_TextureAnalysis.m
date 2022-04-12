% Texture Analysis Toolbox
% Copyright (C) 2014 
% Version 1.0, May/2014 
% Rouzbeh Maani
%__________________________________________________________________________
function spm_TextureAnalysis()
fnlist = {
    'VGLCM TOP 3D',             'spm_texture_glcmtop_module'
    'VGLCM 3D',              'spm_texture_3dglcm_module'
    'GLCM 3D ROI Analysis', 'spm_texture_3dglcm_roi_module'
    'Mask in Original Space', 'spm_mask_module'
};

str = sprintf('%s|', fnlist{:, 1});
str = str(1:(end-1));

fun = spm_input('Texture Analysis',1,'m', str, strvcat(fnlist(:, 2)));
  
eval(fun);

end