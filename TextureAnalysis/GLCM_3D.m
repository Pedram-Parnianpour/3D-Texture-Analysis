function GLCM_3D(varargin)
% GLCM on 3D space
[fileList,maskFiles,OPT,SPACE] = ParseParams(varargin);
tlboxDir = fileparts(mfilename('fullpath'));
currentDir = pwd;
MNIMaskVol = load_untouch_nii([tlboxDir '/data/T1_Mask.nii']);
maskMNI = MNIMaskVol.img>0;
fprintf('Starting Time: %s\n',datestr(now));
for i=1:length(fileList)
    if(SPACE == 0)
        space = 'MNI';
        maskFile = maskFiles{1};
    else
        space = 'Original';
        maskFile = maskFiles{i};
    end    
    MaskVol = load_untouch_nii(maskFile);
    mask = MaskVol.img>0;   
    vol = load_untouch_nii(fileList{i});
    vol = single(vol.img);
    if(size(vol)~= size(mask))
        errordlg(['The size of the mask (' maskFile ') and the volume (' fileList{i} ') is not the same.!']);
        finish;    
    elseif(SPACE == 0 && (size(maskMNI,1)~= size(mask,1)|size(maskMNI,2)~= size(mask,2)|size(maskMNI,3)~= size(mask,3)))
        errordlg(['The size of the mask (' maskFile ') looks different from size of MNI mask[' num2str(size(maskMNI)) '] and the MNI space is used for computation! Please choose a right data or choose the space of computation appropriately.']);
        finish;            
    end
    [d1,d2,d3] = ind2sub(size(mask), find(mask));
    InsideRange = [min(d1) max(d1);min(d2) max(d2);min(d3) max(d3)];
    mask = mask(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
    vol = vol(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
    timerStart = tic;
    fprintf('================\n');
    fprintf('Subject %d: %s \n',i,fileList{i});
    % Intensity Normalization: Assume that the higher/lower 1% are outliers for normalization
    NoisePercent = 0.01;
    I = vol(mask(:));
    sortI = sort(I);
    norm_fact_U_Ind = round((1-NoisePercent)*length(I));
    norm_fact_L_Ind = round(NoisePercent*length(I));
    norm_fact_U = sortI(norm_fact_U_Ind);
    norm_fact_L = sortI(norm_fact_L_Ind);
    newLower = (1/(OPT.quantLevel-1));
    newUpper = 1-newLower;
    I = (I-norm_fact_L)/(norm_fact_U-norm_fact_L)*(newUpper-newLower)+newLower;
    I(I>1) = 1;
    I(I<0) = 0;
    vol(mask) = I;
    vol(~mask(:)) = 0;   
    % End of Intensity Normalization
    [outpath outName] = fileparts(fileList{i});
    glcmFile = [outpath '\' outName '_3D_D' num2str(OPT.D) '_N' num2str(OPT.NeighborSize) '_Q' num2str(OPT.quantLevel) '_' space '.mat']; 
    if (exist(glcmFile, 'file') == 2)
        GLCMS = load(glcmFile);
        GLCMS = GLCMS.GLCMS;
    else
        GLCMS = CreateGLCM_Local(vol,OPT.quantLevel,[0 1],OPT.D,OPT.NeighborSize,mask,SPACE); 
        save(glcmFile,'GLCMS');
    end
    texture = computeGLCMLocalFeat(GLCMS,mask,OPT.glcm_properties);        
    time_total = toc(timerStart);
    fprintf('Total time: %2.2f\n',time_total);
    if(OPT.SmoothKerSize~=0)
        h = fspecial3('gaussian',OPT.SmoothKerSize);
    end    
    vol = load_untouch_nii(fileList{i});
    [outpath outName] = fileparts(fileList{i});
    NormFacts = getNormFact(OPT.glcm_properties,OPT.quantLevel);
    if(SPACE == 1)
        fprintf('Transform Texture Features into MNI\n');
    end        
    for prop=1:length(OPT.glcm_properties)
        img = squeeze(texture(:,:,:,prop));
        if(OPT.SmoothKerSize~=0)
            img = imfilter(img,h);
        end        
        img = (img-min(img(:)))*NormFacts(prop);
        img(~mask) = 0;
        vol.img = zeros(size(vol.img));
        vol.img(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2)) = img;
        vol.img(MaskVol.img<=0) = 0;
        vol.hdr.dime.glmax = max(img(:));
        vol.hdr.dime.glmin = min(img(:));
        vol.hdr.hist.descrip = 'Texture';
        TextureFileName  = [outName '_' OPT.glcm_properties{prop} '_3D_D' num2str(OPT.D) '_N' num2str(OPT.NeighborSize) '_Q' num2str(OPT.quantLevel) '_S' num2str(OPT.SmoothKerSize) '_' space '.nii'];
        save_untouch_nii(vol,[outpath '/' TextureFileName]);  
        if(SPACE == 1)
        % Convert to MNI space
            DeformMap = [outpath '/y_r' outName '.nii'];
            batchName = [tlboxDir '/data/batch_Mask.mat'];
            batchFile = batchName;
            BATCH = load(batchFile);
            BATCH.matlabbatch{1}.spm.util.defs.ofname = outName;
            BATCH.matlabbatch{1}.spm.util.defs.comp{1}.def{1} = DeformMap;
            oldField = BATCH.matlabbatch{1}.spm.util.defs.fnames;
            oldField = strrep(oldField, 'T1_Mask.nii', TextureFileName);
            newField = strrep(oldField, 'D:\Academic\Research\MyCode\VoxelwiseTextureAnalysis\data', outpath);
            BATCH.matlabbatch{1}.spm.util.defs.fnames = newField;
            spm_jobman('run',BATCH.matlabbatch);
            movefile([currentDir '/w' TextureFileName],[outpath '/' TextureFileName],'f');
            delete([currentDir '/y_' outName '.nii']);
            volToSave = load_untouch_nii([outpath '/' TextureFileName]);
            volToSave.img(~maskMNI)=0;
            save_untouch_nii(volToSave,[outpath '/' TextureFileName]); 
        end        
    end
end
fprintf('End Time: %s\n',datestr(now));
end
%--------------------------------------------------------------------------
function offset = ComputeOffsets(dis)
    [delta1,delta2,delta3] = ndgrid(-dis:dis,-dis:dis,-dis:dis);
    inside = round(sqrt(delta1.*delta1+delta2.*delta2+delta3.*delta3))<=dis;
    offset = [delta1(inside(:)),delta2(inside(:)),delta3(inside(:))];
    offset(offset(:,1)==0&offset(:,2)==0&offset(:,3)==0,:) = [];
    % For symmetric glcm we need half of offsets
    % offset = offset(1:size(offset,1)/2,:);
end
%--------------------------------------------------------------------------


function texture = computeGLCMLocalFeat(glcm,mask,GLCM_feat)
texture = zeros(size(glcm,1),size(glcm,2),size(glcm,3),length(GLCM_feat),'single');
fprintf('Computing Features: ')
n=0;
for i=1:size(glcm,1)
    for j=1:size(glcm,2)
        for k=1:size(glcm,3)
            if(mask(i,j,k)==0)
                continue;
            end
            glcm_norm = squeeze(glcm(i,j,k,:,:));            
            st = computeFeature(glcm_norm,GLCM_feat); 
            texture(i,j,k,:) = st;
        end
    end
    fprintf(repmat('\b',1,n));
    msg = sprintf('%2.2f',i/size(glcm,1)*100);
    fprintf([msg '%%']);
    n=numel(msg)+1;                    
end
fprintf('\n');
end

function feature = computeFeature(glcm,GLCM_feat_all)
feature = zeros(length(GLCM_feat_all),1);
size_glcm_1 = size(glcm,1); 
size_glcm_2 = size(glcm,2); 
Pij = glcm;
[i,j] = meshgrid(1:size_glcm_1,1:size_glcm_2); 
p_x = squeeze(sum(Pij,2)); 
p_y = squeeze(sum(Pij,1))'; 
glcm_mean = mean(Pij(:)); 
idx1 = (i+j)-1; 
p_xplusy = zeros((2*size_glcm_1 - 1),1);
for aux = 1:max(idx1(:)) 
    p_xplusy(aux) = sum(Pij(idx1==aux)); 
end 
ii = (1:(2*size_glcm_1-1))'; 
jj = (0:size_glcm_1-1)';
idx2 = abs(i-j)+1; 
p_xminusy = zeros((size_glcm_1),1);
for aux = 1:max(idx2(:)) 
    p_xminusy(aux) = sum(Pij(idx2==aux)); 
end 
u_x = sum(sum(i.*Pij)); 
u_y = sum(sum(j.*Pij));     
s_x = sum(sum(Pij.*((i-u_x).^2)))^0.5; 
s_y = sum(sum(Pij.*((j-u_y).^2)))^0.5; 

for prop=1:length(GLCM_feat_all)
    GLCM_feat = GLCM_feat_all{prop};
    % Autocorrelation
    if(strcmp(GLCM_feat,'autoc')==1)    
        feature(prop) = sum(sum(Pij.*(i.*j))); 
    % Contrast 
    elseif(strcmp(GLCM_feat,'contr')==1)
        feature(prop) = sum(sum((abs(i-j).^2).*Pij)); 
    % Dissimilarity 
    elseif(strcmp(GLCM_feat,'dissi')==1)
        feature(prop) = sum(sum(abs(i-j).*Pij)); 
    % Energy 
    elseif(strcmp(GLCM_feat,'energ')==1)
        feature(prop) = sum(sum(Pij.^2)); 
    % Entropy 
    elseif(strcmp(GLCM_feat,'entro')==1)
        feature(prop) = -sum(sum(Pij.*log(Pij+eps))); 
    % Homogeneity Matlab 
    elseif(strcmp(GLCM_feat,'homom')==1)
        feature(prop) =  sum(sum(Pij./(1+abs(i-j)))); 
    % Homogeneity Paper 
    elseif(strcmp(GLCM_feat,'homop')==1)
        feature(prop) = sum(sum(Pij./(1+abs(i-j).^2))); 
    % Sum of squares: Variance 
    elseif(strcmp(GLCM_feat,'sosvh')==1)
        feature(prop) = sum(sum(Pij.*((j-glcm_mean).^2))); 
    % Inverse difference normalized 
    elseif(strcmp(GLCM_feat,'indnc')==1)
        feature(prop) = sum(sum(Pij./(1+(abs(i-j)./size_glcm_1)))); 
    % Inverse difference moment normalized 
    elseif(strcmp(GLCM_feat,'idmnc')==1)
        feature(prop) = sum(sum(Pij./(1+((i-j)./size_glcm_1).^2))); 
    % Maximum probability 
    elseif(strcmp(GLCM_feat,'maxpr')==1)
        feature(prop) = max(Pij(:)); 
    % Sum average 
    elseif(strcmp(GLCM_feat,'savgh')==1)
        feature(prop) = sum((ii+1).*p_xplusy); 
    % Sum entropy 
    elseif(strcmp(GLCM_feat,'senth')==1)
        feature(prop) = -sum(p_xplusy.*log(p_xplusy+eps));     
    % Sum variance 
    elseif(strcmp(GLCM_feat,'svarh')==1)
        senth = -sum(p_xplusy.*log(p_xplusy+eps)); 
        feature(prop) = sum((((ii+1) - senth).^2).*p_xplusy); 
    % Difference entropy 
    elseif(strcmp(GLCM_feat,'denth')==1)
        feature(prop) = -sum(p_xminusy.*log(p_xminusy+eps)); 
    % Difference variance 
    elseif(strcmp(GLCM_feat,'dvarh')==1)
        feature(prop) = sum((jj.^2).*p_xminusy); 
    % Correlation Matlab     
    elseif(strcmp(GLCM_feat,'corrm')==1)
        corm = sum(sum(Pij.*(i-u_x).*(j-u_y))); 
        feature(prop) = corm/(s_x*s_y); 
    % Correlation paper 
    elseif(strcmp(GLCM_feat,'corrp')==1)
        corp = sum(sum(Pij.*(i.*j))); 
        feature(prop) = (corp-u_x*u_y)/(s_x*s_y); 
    % Cluster Prominence 
    elseif(strcmp(GLCM_feat,'cprom')==1)
        feature(prop) = sum(sum(Pij.*((i+j-u_x-u_y).^4))); 
    % Cluster Shade    
    elseif(strcmp(GLCM_feat,'cshad')==1)
        feature(prop) = sum(sum(Pij.*((i+j-u_x-u_y).^3)));   
    % Information measure of correlation 1 
    elseif(strcmp(GLCM_feat,'inf1h')==1)
        hx = -sum(p_x.*log(p_x+eps)); 
        hy = -sum(p_y.*log(p_y+eps)); 
        hxy = -sum(sum(Pij.*log(Pij+eps))); 
        hxy1 = -sum(sum(Pij.*log(p_x*p_y' + eps))); 
        feature(prop) = (hxy-hxy1)/(max([hx,hy])); 
    % Information measure of correlation 2
    elseif(strcmp(GLCM_feat,'inf2h')==1)
        hxy = -sum(sum(Pij.*log(Pij+eps))); 
        hxy2 = -sum(sum((p_x*p_y').*log(p_x*p_y' + eps))); 
        feature(prop) = (1-exp(-2*(hxy2-hxy)))^0.5;     
    end  
end
feature(isnan(feature)) = 0;
feature(isinf(feature)) = 0;
end

%--------------------------------------------------------------------------
function [GLCMS] = CreateGLCM_Local(I, NL, GL,D,NeighborSizeRad,mask,SPACE)
if GL(2) == GL(1)
    SI = ones(size(I));
else
    slope = (NL-1) / (GL(2) - GL(1));
    intercept = 1 - (slope*(GL(1)));
    SI = round(imlincomb(slope,I,intercept,'double'));
end
SI(SI > NL) = NL;
SI(SI < 1) = 1;
vol_Gray = single(SI);
SIZ = size(I);
clear I SI;
if(SPACE==0)
    mappingFile = ['Mapping_D' num2str(D) '_N' num2str(NeighborSizeRad) '.mat'];
    if (exist(mappingFile, 'file') == 2)
        Mapping = load(mappingFile);
        OffsetsMatrixS = Mapping.OffsetsMatrixS;
        OffsetsMatrixE = Mapping.OffsetsMatrixE;
    else
        [OffsetsMatrixS OffsetsMatrixE] = AllOffsetsAllNeighbors(vol_Gray,D,NeighborSizeRad,mask);
        save(mappingFile,'OffsetsMatrixS','OffsetsMatrixE','-v7.3');
    end
else
    [OffsetsMatrixS OffsetsMatrixE] = AllOffsetsAllNeighbors(vol_Gray,D,NeighborSizeRad,mask);
end
InsideMask = OffsetsMatrixS>0&OffsetsMatrixE>0;
glcm = zeros(length(InsideMask),NL,NL,'single');
fprintf('Computing GLCM: ')
n=0;
for i=1:NL
    for j=1:NL
        G = zeros(size(InsideMask),'uint8');
        G(InsideMask) = vol_Gray(OffsetsMatrixS(InsideMask))==i&vol_Gray(OffsetsMatrixE(InsideMask))==j;
        glcm(:,i,j) = sum(sum(G,3),2);
    end
    fprintf(repmat('\b',1,n));
    msg = sprintf('%2.2f',i/NL*100);
    fprintf([msg '%%']);
    n=numel(msg)+1;                    
end
fprintf('\n');
clear vol_Gray;
norm_fact = sum(sum(glcm,3),2);
glcm = glcm./repmat(norm_fact,[1 NL NL]);
mask = repmat(mask,[1 1 1 NL NL]);
GLCMS = zeros(SIZ(1),SIZ(2),SIZ(3),NL,NL,'single');
GLCMS(mask) = glcm;
end

function [OffsetsMatrixS OffsetsMatrixE] = AllOffsetsAllNeighbors(I,D,NeighborSizeRad,mask)
Y = single(1:size(I,1));
X = single(1:size(I,2));
Z = single(1:size(I,3));
[X,Y,Z] = meshgrid(X,Y,Z);
X = X(mask(:));
Y = Y(mask(:));
Z = Z(mask(:));
IND_s = single(sub2ind(size(I),Y,X,Z));
% Computing all offsets at the voxel
Offsets = single(ComputeOffsets(D));
% For speed we can use half of offsets
% Offsets = Offsets(1:size(Offsets,1)/2,:);
Neighbors = single(ComputeOffsets(NeighborSizeRad));
OffsetsMatrixS = single(repmat(IND_s,[1 ,size(Offsets,1),size(Neighbors,1)+1]));
OffsetsMatrixE = zeros(size(OffsetsMatrixS),'single');
fprintf('Computing Maps: ');
n=0;
for ofst=1:size(Offsets,1)
    offset = Offsets(ofst,:);
    Y_e = Y+offset(1);
    X_e = X+offset(2);
    Z_e = Z+offset(3);
    isInMask = Y_e>0&X_e>0&Z_e>0&Y_e<=size(I,1)&X_e<=size(I,2)&Z_e<=size(I,3);
    IND_e = sub2ind(size(I),Y_e(isInMask),X_e(isInMask),Z_e(isInMask));%End point of vector
    %Check weather it is masked
    isInMask2 = mask(IND_e);
    isInMask(isInMask) = isInMask2;
    OffsetsMatrixE(isInMask,ofst,1) = IND_e(isInMask2);
    % Computing all offsets at the neighbors of the voxel
    for nb=1:size(Neighbors,1)
        NeighborOffset = Neighbors(nb,:);
        Y_ns = Y+NeighborOffset(1);
        X_ns = X+NeighborOffset(2);
        Z_ns = Z+NeighborOffset(3);
        isInMask = Y_ns>0&X_ns>0&Z_ns>0&Y_ns<=size(I,1)&X_ns<=size(I,2)&Z_ns<=size(I,3);
        IND_ns = sub2ind(size(I),Y_ns(isInMask),X_ns(isInMask),Z_ns(isInMask));
        isInMask2 = mask(IND_ns);
        isInMask(isInMask) = isInMask2;
        OffsetsMatrixS(isInMask,ofst,nb+1) = IND_ns(isInMask2);
        %Check weather it is masked
        Y_ne = Y_e+NeighborOffset(1);
        X_ne = X_e+NeighborOffset(2);
        Z_ne = Z_e+NeighborOffset(3);
        isInMask = Y_ne>0&X_ne>0&Z_ne>0&Y_ne<=size(I,1)&X_ne<=size(I,2)&Z_ne<=size(I,3);
        IND_ne = sub2ind(size(I),Y_ne(isInMask),X_ne(isInMask),Z_ne(isInMask));
        isInMask2 = mask(IND_ne);
        isInMask(isInMask) = isInMask2;
        OffsetsMatrixE(isInMask,ofst,nb+1) = IND_ne(isInMask2);
    end
    fprintf(repmat('\b',1,n));
    msg = sprintf('%2.2f',ofst/size(Offsets,1)*100);
    fprintf([msg '%%']);
    n=numel(msg)+1;                
end
fprintf('\n');
end
%--------------------------------------------------------------------------
function NormFacts = getNormFact(GLCM_feat_all,Q)
NormFacts = zeros(length(GLCM_feat_all),1);
for prop=1:length(GLCM_feat_all)
    GLCM_feat = GLCM_feat_all{prop};
    % Autocorrelation
    if(strcmp(GLCM_feat,'autoc')==1)    
        NormFacts(prop) = 16/Q; 
    % Contrast 
    elseif(strcmp(GLCM_feat,'contr')==1)
        NormFacts(prop) = 5^((16/Q)-1)*4; 
    % Dissimilarity 
    elseif(strcmp(GLCM_feat,'dissi')==1)
        NormFacts(prop) = 50; 
    % Energy 
    elseif(strcmp(GLCM_feat,'energ')==1)
        NormFacts(prop) = 200; 
    % Entropy 
    elseif(strcmp(GLCM_feat,'entro')==1)
        NormFacts(prop) = 50; 
    % Homogeneity Matlab 
    elseif(strcmp(GLCM_feat,'homom')==1)
        NormFacts(prop) =  200; 
    % Homogeneity Paper 
    elseif(strcmp(GLCM_feat,'homop')==1)
        NormFacts(prop) = 200; 
    % Sum of squares: Variance 
    elseif(strcmp(GLCM_feat,'sosvh')==1)
        NormFacts(prop) = 4^((16/Q)-1); 
    % Inverse difference normalized 
    elseif(strcmp(GLCM_feat,'indnc')==1)
        NormFacts(prop) = 200; 
    % Inverse difference moment normalized 
    elseif(strcmp(GLCM_feat,'idmnc')==1)
        NormFacts(prop) = 200; 
    % Maximum probability 
    elseif(strcmp(GLCM_feat,'maxpr')==1)
        NormFacts(prop) = 200; 
    % Sum average 
    elseif(strcmp(GLCM_feat,'savgh')==1)
        NormFacts(prop) = 2^(16/Q)*4; 
    % Sum entropy 
    elseif(strcmp(GLCM_feat,'senth')==1)
        NormFacts(prop) = 100;     
    % Sum variance 
    elseif(strcmp(GLCM_feat,'svarh')==1)
        NormFacts(prop) = (0.25)^((Q/8)-1); 
    % Difference entropy 
    elseif(strcmp(GLCM_feat,'denth')==1)
        NormFacts(prop) = 100; 
    % Difference variance 
    elseif(strcmp(GLCM_feat,'dvarh')==1)
        NormFacts(prop) = 5^((16/Q)-1)*4; 
    % Correlation Matlab     
    elseif(strcmp(GLCM_feat,'corrm')==1)
        NormFacts(prop) = 100; 
    % Correlation paper 
    elseif(strcmp(GLCM_feat,'corrp')==1)
        NormFacts(prop) = 1; 
    % Cluster Prominence 
    elseif(strcmp(GLCM_feat,'cprom')==1)
        if(Q<=8)
            NormFacts(prop) = 1;
        else
            NormFacts(prop) = 0.1^(log2((Q/8))+1); 
        end
    % Cluster Shade    
    elseif(strcmp(GLCM_feat,'cshad')==1)
        NormFacts(prop) = 4^((16/Q)-1)/(4^log2(Q/8));   
    % Information measure of correlation 1 
    elseif(strcmp(GLCM_feat,'inf1h')==1)
        NormFacts(prop) = 400; 
    % Information measure of correlation 2
    elseif(strcmp(GLCM_feat,'inf2h')==1)
        NormFacts(prop) = 400;     
    end  
end
end

function [fileList,maskFiles,OPT,SPACE] = ParseParams(varargin)
fileList = cell(1,length(varargin{1}{1}.data));
for i=1:length(fileList)
    fileList{i} = strtok(varargin{1}{1}.data{i},',');
end
maskFiles = cell(1,length(varargin{1}{1}.masks));
for i=1:length(maskFiles)
    maskFiles{i} = strtok(varargin{1}{1}.masks{i},',');
end
OPT.D = varargin{1}{1}.opts.distance;
OPT.NeighborSize = varargin{1}{1}.opts.neighborsize;
OPT.quantLevel = varargin{1}{1}.opts.quantlevel;
OPT.SmoothKerSize = varargin{1}{1}.kernelsize;
SPACE = varargin{1}{1}.space;
if(SPACE==0 && length(maskFiles)~=1)
    errordlg('The number of mask(s) should be 1 for computation in MNI space!');
    finish;
end
if(SPACE==1 && length(maskFiles)~=length(fileList))
    errordlg('The number of mask(s) should be the same as the number of volumes for computation in Original space!');
    finish;
end
glcm_properties = {'autoc','contr','corrm','corrp','cprom','cshad','dissi','energ','entro','homom','homop','maxpr','sosvh','savgh','svarh','senth','dvarh','denth','inf1h','inf2h','indnc','idmnc'};
OPT.glcm_properties = {};
r=1;
for i=1:length(glcm_properties)
    [~,present] = evalc(['varargin{1}{1}.features.' glcm_properties{i}]);
    if(present==1)
        OPT.glcm_properties{r} = glcm_properties{i};
        r = r+1;
    end
end

end
