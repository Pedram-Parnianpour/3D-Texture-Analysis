function TA = tbx_cfg_TextureAnalysis()
% Configuration file for Texture Analysis


addpath(fileparts(which(mfilename)));

%==========================================================================
%============================== GLCM TOP ==================================
%==========================================================================

%------------------------------------------------------------------------
% Input files
%------------------------------------------------------------------------

data = cfg_files;
data.tag  = 'data';
data.name = 'Volumes';
data.help = {'Select raw data (e.g. T1 images) for texture alanysis.'};
data.filter = 'image';
data.ufilter = '.*';
data.num     = [1 Inf];

%------------------------------------------------------------------------
% Inverse Deformation Maps files
%------------------------------------------------------------------------

invdefmap = cfg_files;
invdefmap.tag  = 'invdefmap';
invdefmap.name = 'Inv Deformation Maps';
invdefmap.help = {'Select inverse deformation maps for mask (with prefix iy_r) computation in original space.'};
invdefmap.filter = 'image';
invdefmap.ufilter = 'iy_r';
invdefmap.num     = [1 Inf];

%------------------------------------------------------------------------
% Mask file
%------------------------------------------------------------------------

mask = cfg_files;
mask.tag  = 'mask';
mask.name = 'Mask';
mask.help = {'Select the mask (e.g. mask of intracranial tissues) for texture alanysis.'};
mask.filter = 'image';
mask.ufilter = '.*';
mask.def  =  @(val)cfg_TA_get_defaults('mask.def'); 
mask.num     = [1 1];

%------------------------------------------------------------------------
% Masks file (for ROI analysis when each volume uses a mask)
%------------------------------------------------------------------------

masks = cfg_files;
masks.tag  = 'masks';
masks.name = 'Mask(s)';
masks.help = {'Select the mask or masks for ROI alanysis.',...
              'You must choose either ONE mask to apply to all images ',...
              '(when they are registered into the same template) OR one ',...
              'mask for each volume. The order of the masks should be the same ',...
              'as the order of volumes.'};
masks.filter = 'image';
masks.ufilter = '.*';
masks.num     = [1 Inf];

% ---------------------------------------------------------------------
% dir Directory (for writing ROI excel file)
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {'Select a directory where the excel file containing the texture values of ROI analysis will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

%------------------------------------------------------------------------
% Parameters for texture analysis
%------------------------------------------------------------------------

quantlevel    = cfg_menu;
quantlevel.tag = 'quantlevel';
quantlevel.name = 'Quantize Level';
quantlevel.labels = {'8','16'};
quantlevel.values = {8,16};
quantlevel.def  =  @(val)cfg_TA_get_defaults('quantlevel.def');
quantlevel.help = {[...
'The number of gray levels to compute GLCM. ',...
'In most of the cases 8 is enough. Increasing the level to 16 may improve the results ',...
'with the expense of more processing time.']};

%------------------------------------------------------------------------

space    = cfg_menu;
space.tag = 'space';
space.name = 'Computation Space';
space.labels = {'MNI','Original'};
space.values = {0,1};
space.def  =  @(val)cfg_TA_get_defaults('space.def');
space.help = {[...
'The space in which texture features are computed. ',...
'The space could be the original space (before registration) or the MNI space (after registration).']};

%------------------------------------------------------------------------

smoothing    = cfg_menu;
smoothing.tag = 'kernelsize';
smoothing.name = 'Smoothing Kernel Size';
smoothing.labels = {'0','3','5'};
smoothing.values = {0,3,5};
smoothing.def  =  @(val)cfg_TA_get_defaults('smoothing.def');
smoothing.help = {[...
'The size of the smoothing Guassian kernel to apply to the texture features. ',...
'The default value is 0 because in general the texture are smooth. To further smooth ',...
'and remove unexpected noise use this option.']};

%------------------------------------------------------------------------

neighborsize    = cfg_menu;
neighborsize.tag = 'neighborsize';
neighborsize.name = 'Neighborhood Radius';
neighborsize.labels = {'1','2','3'};
neighborsize.values = {1,2,3};
neighborsize.def  = @(val)cfg_TA_get_defaults('neighborsize.def');
neighborsize.help = {[...
'The radius of the neighborhood around each voxel to compute GLCM. ',...
'Depending on the scale of pathological changes one may choose the radius ',...
'(e.g., for changes occuring at small scales choose 1 voxels). ']};

%------------------------------------------------------------------------

distance = cfg_menu;
distance.tag = 'distance';
distance.name = 'Distance';
distance.labels = {'1','2','3','4'};
distance.values = {1,2,3,4};
distance.def  = @(val)cfg_TA_get_defaults('distance.def');
distance.help = {[...
'The largest distance to compute GLCM. The the offsets smaller than or equal ',...
'to this distance will be considered to compute GLCM.']};

%------------------------------------------------------------------------

opts      = cfg_branch;
opts.tag = 'opts';
opts.name = 'VGLCM Parameters';
opts.val = {quantlevel,neighborsize,distance};
opts.help = {[...
'The options can be adjusted in order to improve the performance of texture ',...
'analysis.  Choosing the best parameters is a matter ',...
'of empirical experiments. For instance for changes occuring at small scales the ',...
'distance and the neighboring size should be considered small. The default values ',...
'are good points to start.']};

%------------------------------------------------------------------------

optsROI      = cfg_branch;
optsROI.tag = 'optsROI';
optsROI.name = 'GLCM Parameters';
optsROI.val = {quantlevel,distance};
optsROI.help = {[...
'The options can be adjusted in order to improve the performance of texture ',...
'analysis.  Choosing the best parameters depends of the data. ',...
'For instance for changes occuring at small scales the ',...
'neighboring diameter should be considered small. The default values ',...
'are good points to start.']};
%-----------------------------------------------------------------------
% Texture Features
%-----------------------------------------------------------------------

autoc    = cfg_menu;
autoc.tag = 'autoc';
autoc.name = 'Autocorrelation';
autoc.labels = {'Yes','No'};
autoc.values = {1,0};
autoc.def  = @(val)cfg_TA_get_defaults('autoc.def');
autoc.help = {'Autocorrelation feature of GLCM.'};

contr    = cfg_menu;
contr.tag = 'contr';
contr.name = 'Contrast';
contr.labels = {'Yes','No'};
contr.values = {1,0};
contr.def  = @(val)cfg_TA_get_defaults('contr.def');
contr.help = {'Contrast feature of GLCM.'};

dissi    = cfg_menu;
dissi.tag = 'dissi';
dissi.name = 'Dissimilarity';
dissi.labels = {'Yes','No'};
dissi.values = {1,0};
dissi.def  = @(val)cfg_TA_get_defaults('dissi.def');
dissi.help = {'Dissimilarity feature of GLCM.'};

energ    = cfg_menu;
energ.tag = 'energ';
energ.name = 'Energy';
energ.labels = {'Yes','No'};
energ.values = {1,0};
energ.def  = @(val)cfg_TA_get_defaults('energ.def');
energ.help = {'Energy feature of GLCM.'};

entro    = cfg_menu;
entro.tag = 'entro';
entro.name = 'Entropy';
entro.labels = {'Yes','No'};
entro.values = {1,0};
entro.def  = @(val)cfg_TA_get_defaults('entro.def');
entro.help = {'Entropy feature of GLCM.'};

homom    = cfg_menu;
homom.tag = 'homom';
homom.name = 'Homogeneity(1)';
homom.labels = {'Yes','No'};
homom.values = {1,0};
homom.def  = @(val)cfg_TA_get_defaults('homom.def');
homom.help = {'Homogeneity feature of GLCM using the formula in Matlab.'};

homop    = cfg_menu;
homop.tag = 'homop';
homop.name = 'Homogeneity(2)';
homop.labels = {'Yes','No'};
homop.values = {1,0};
homop.def  = @(val)cfg_TA_get_defaults('homop.def');
homop.help = {'Homogeneity feature of GLCM using the formula in the GLCM paper.'};

sosvh    = cfg_menu;
sosvh.tag = 'sosvh';
sosvh.name = 'Sum of squares: Variance';
sosvh.labels = {'Yes','No'};
sosvh.values = {1,0};
sosvh.def  = @(val)cfg_TA_get_defaults('sosvh.def');
sosvh.help = {'Sum of squares: Variance feature of GLCM.'};

indnc    = cfg_menu;
indnc.tag = 'indnc';
indnc.name = 'Inverse difference normalized';
indnc.labels = {'Yes','No'};
indnc.values = {1,0};
indnc.def  = @(val)cfg_TA_get_defaults('indnc.def');
indnc.help = {'Inverse difference normalized feature of GLCM.'};

idmnc    = cfg_menu;
idmnc.tag = 'idmnc';
idmnc.name = 'Inverse difference moment normalized';
idmnc.labels = {'Yes','No'};
idmnc.values = {1,0};
idmnc.def  = @(val)cfg_TA_get_defaults('idmnc.def');
idmnc.help = {'Inverse difference moment normalized feature of GLCM.'};

maxpr    = cfg_menu;
maxpr.tag = 'maxpr';
maxpr.name = 'Maximum probability';
maxpr.labels = {'Yes','No'};
maxpr.values = {1,0};
maxpr.def  = @(val)cfg_TA_get_defaults('maxpr.def');
maxpr.help = {'Maximum probability feature of GLCM.'};

savgh    = cfg_menu;
savgh.tag = 'savgh';
savgh.name = 'Sum average';
savgh.labels = {'Yes','No'};
savgh.values = {1,0};
savgh.def  = @(val)cfg_TA_get_defaults('savgh.def');
savgh.help = {'Sum average feature of GLCM.'};

senth    = cfg_menu;
senth.tag = 'senth';
senth.name = 'Sum entropy';
senth.labels = {'Yes','No'};
senth.values = {1,0};
senth.def  = @(val)cfg_TA_get_defaults('senth.def');
senth.help = {'Sum entropy feature of GLCM.'};

svarh    = cfg_menu;
svarh.tag = 'svarh';
svarh.name = 'Sum variance';
svarh.labels = {'Yes','No'};
svarh.values = {1,0};
svarh.def  = @(val)cfg_TA_get_defaults('svarh.def');
svarh.help = {'Sum variance feature of GLCM.'};

denth    = cfg_menu;
denth.tag = 'denth';
denth.name = 'Difference entropy';
denth.labels = {'Yes','No'};
denth.values = {1,0};
denth.def  = @(val)cfg_TA_get_defaults('denth.def');
denth.help = {'Difference entropy feature of GLCM.'};

dvarh    = cfg_menu;
dvarh.tag = 'dvarh';
dvarh.name = 'Difference variance';
dvarh.labels = {'Yes','No'};
dvarh.values = {1,0};
dvarh.def  = @(val)cfg_TA_get_defaults('dvarh.def');
dvarh.help = {'Difference variance feature of GLCM.'};

corrm    = cfg_menu;
corrm.tag = 'corrm';
corrm.name = 'Correlation(1)';
corrm.labels = {'Yes','No'};
corrm.values = {1,0};
corrm.def  = @(val)cfg_TA_get_defaults('corrm.def');
corrm.help = {'Correlation feature of GLCM using the formula in Matlab.'};

corrp    = cfg_menu;
corrp.tag = 'corrp';
corrp.name = 'Correlation(2)';
corrp.labels = {'Yes','No'};
corrp.values = {1,0};
corrp.def  = @(val)cfg_TA_get_defaults('corrp.def');
corrp.help = {'Correlation feature of GLCM using the formula in the GLCM paper.'};

cprom    = cfg_menu;
cprom.tag = 'cprom';
cprom.name = 'Cluster Prominence';
cprom.labels = {'Yes','No'};
cprom.values = {1,0};
cprom.def  = @(val)cfg_TA_get_defaults('cprom.def');
cprom.help = {'Cluster Prominence feature of GLCM.'};

cshad    = cfg_menu;
cshad.tag = 'cshad';
cshad.name = 'Cluster Shade';
cshad.labels = {'Yes','No'};
cshad.values = {1,0};
cshad.def  = @(val)cfg_TA_get_defaults('cshad.def');
cshad.help = {'Cluster Shade feature of GLCM.'};

inf1h    = cfg_menu;
inf1h.tag = 'inf1h';
inf1h.name = 'Information measure of correlation 1';
inf1h.labels = {'Yes','No'};
inf1h.values = {1,0};
inf1h.def  = @(val)cfg_TA_get_defaults('inf1h.def');
inf1h.help = {'Information measure of correlation 1 feature of GLCM.'};

inf2h    = cfg_menu;
inf2h.tag = 'inf2h';
inf2h.name = 'Information measure of correlation 2';
inf2h.labels = {'Yes','No'};
inf2h.values = {1,0};
inf2h.def  = @(val)cfg_TA_get_defaults('inf2h.def');
inf2h.help = {'Information measure of correlation 2 feature of GLCM.'};

features      = cfg_branch;
features.tag = 'features';
features.name = 'Texture Features';
features.val = {autoc,contr,corrm,corrp,cprom,cshad,dissi,denth,dvarh,energ,entro,homom,homop,inf1h,inf2h,indnc,idmnc,maxpr,savgh,senth,sosvh,svarh};
features.help = {'Choose the texture features to analysis.'};

%--------------------------------------------------------------------------
% Module Definition
%--------------------------------------------------------------------------
%==========================================================================
%============================= VGLCM TOP 3D ===============================
%==========================================================================

glcmtop      = cfg_exbranch;
glcmtop.tag = 'glcmtop';
glcmtop.name = 'Texture Analysis: VGLCM TOP 3D';
glcmtop.val = {data,masks,opts,space,smoothing,features};
glcmtop.prog   = @GLCM_TOP;
glcmtop.help   = {'Voxel-based GLCM on Three Orthogonal Planes in 3D space (VGLCM-TOP-3D).'};

%==========================================================================
%============================== VGLCM 3D ==================================
%==========================================================================

glcm3d      = cfg_exbranch;
glcm3d.tag = 'glcm3d';
glcm3d.name = 'Texture Analysis: VGLCM 3D';
glcm3d.val = {data,masks,opts,space,smoothing,features};
glcm3d.prog   = @GLCM_3D;
glcm3d.help   = {'Voxel-based GLCM in 3D space!'};


%==========================================================================
%============================= 3D GLCM ROI ================================
%==========================================================================

glcm3droi      = cfg_exbranch;
glcm3droi.tag = 'glcm3droi';
glcm3droi.name = 'Texture Analysis: GLCM 3D ROI';
glcm3droi.val = {dir,data,masks,optsROI,features};
glcm3droi.prog   = @GLCM_ROI;
glcm3droi.help   = {'Texture Analysis on ROI using 3D GLCM!'};

%==========================================================================
%======================== Mask in Original Space ==========================
%==========================================================================

computemask      = cfg_exbranch;
computemask.tag = 'computemask';
computemask.name = 'Compute Mask of the Original Space';
computemask.val = {invdefmap};
computemask.prog   = @Compute_Mask_OrigSpace;
computemask.help   = {'Computes Mask of the Original Space by Mapping MNI Mask to the Original Space!'};


TA  = cfg_choice;
TA.name = 'Texture Analysis';
TA.tag  = 'TA';
TA.values = {glcmtop,glcm3d,glcm3droi,computemask};


end
