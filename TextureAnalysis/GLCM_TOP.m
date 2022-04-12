function GLCM_TOP(varargin)
% GLCM on Three Orthogonal Planes (TOP)
[fileList,maskFiles,OPT,SPACE] = ParseParams(varargin);
tlboxDir = fileparts(mfilename('fullpath'));
currentDir = pwd;
MNIMaskVol = load_untouch_nii([tlboxDir '/data/T1_Mask.nii']);
maskMNI = MNIMaskVol.img>0;
% MaskVol = load_untouch_nii(maskFile);
% mask = MaskVol.img>0;
% [d1,d2,d3] = ind2sub(size(mask), find(mask));
% InsideRange = [min(d1) max(d1);min(d2) max(d2);min(d3) max(d3)];
% mask = mask(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
offset = ComputeOffsets(OPT.D);
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
    end
    [d1,d2,d3] = ind2sub(size(mask), find(mask));
    InsideRange = [min(d1) max(d1);min(d2) max(d2);min(d3) max(d3)];
    mask = mask(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
    vol = vol(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
    texture = zeros(size(vol,1),size(vol,2),size(vol,3),length(OPT.glcm_properties),3,'single');
    timerStart = tic;
    fprintf('================\n');
    fprintf('Subject %d: %s \n',i,fileList{i});
    fprintf('Texture Features in %s space\n',space);
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
    for v = 1:3%view
        fprintf('Plane %d: ',v);
        n=0;
        for s=1:size(vol,v) 
            fprintf(repmat('\b',1,n));
            msg = [num2str(s) '/' num2str(size(vol,v))];
            fprintf(msg);
            n=numel(msg);            
            [I maskI] = readImgFromVol(vol,mask,v,s);
            I(~maskI) = NaN;
            [d1,d2] = ind2sub(size(maskI), find(maskI));
            InsideRangeI = [min(d1) max(d1);min(d2) max(d2)];
            if(isempty(InsideRangeI))
                continue;
            end
            GLCM = CreateGLCM_Local(I, offset, OPT.quantLevel, [0 1],OPT.NeighborSize,InsideRangeI);
            GLCM = mean(GLCM,5);% Averaging all offsets
            switch(v)
                case 1
                    texture(s,:,:,:,v) = computeGLCMLocalFeat(GLCM,InsideRangeI,OPT.glcm_properties);
                case 2
                    texture(:,s,:,:,v) = computeGLCMLocalFeat(GLCM,InsideRangeI,OPT.glcm_properties);
                case 3
                    texture(:,:,s,:,v) = computeGLCMLocalFeat(GLCM,InsideRangeI,OPT.glcm_properties);
            end
        end
        fprintf('\n');
    end
    % Averaging over views
    texture = mean(texture,5);    
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
        TextureFileName  = [outName '_' OPT.glcm_properties{prop} '_TOP3D_D' num2str(OPT.D) '_N' num2str(OPT.NeighborSize) '_Q' num2str(OPT.quantLevel) '_S' num2str(OPT.SmoothKerSize) '_' space '.nii'];
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
    [delta1,delta2] = ndgrid(-dis:dis,-dis:dis);
    inside = round(sqrt(delta1.*delta1+delta2.*delta2))<=dis;
    offset = [delta1(inside(:)),delta2(inside(:))];
    offset(offset(:,1)==0&offset(:,2)==0,:) = [];
    % For symmetric glcm we need half of offsets
%     offset = offset(1:size(offset,1)/2,:);
end
%--------------------------------------------------------------------------


function texture = computeGLCMLocalFeat(glcm,InsideRange,GLCM_feat)
texture = zeros(size(glcm,1),size(glcm,2),length(GLCM_feat),'single');
glcm = glcm(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),:,:);
for i=1:size(glcm,1)
    for j=1:size(glcm,2)
            g = squeeze(glcm(i,j,:,:));
            norm_fact = repmat(sum(sum(g,1),2),[size(g,1),size(g,2),1]);
            if(sum(norm_fact)==0)
                continue;
            end
            glcm_norm = g./norm_fact;
            glcm_norm(isnan(glcm_norm))=0;
            st = computeFeature(glcm_norm,GLCM_feat); 
            texture(i+InsideRange(1,1)-1,j+InsideRange(2,1)-1,:) = st;
    end
end
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

function [I maskI]=readImgFromVol(vol,mask,view,sliceNO)
switch(view)
    case 1
        I = squeeze(vol(sliceNO,:,:));
        maskI = squeeze(mask(sliceNO,:,:));
    case 2
        I = squeeze(vol(:,sliceNO,:));
        maskI = squeeze(mask(:,sliceNO,:));
    case 3
        I = squeeze(vol(:,:,sliceNO));
        maskI = squeeze(mask(:,:,sliceNO));
end
end
%--------------------------------------------------------------------------
function [GLCMS] = CreateGLCM_Local(I, Offset, NL, GL,NeighborSize,InsideRange)
Iorig = I;
I = I(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2));
if GL(2) == GL(1)
    SI = ones(size(I));
else
    slope = (NL-1) / (GL(2) - GL(1));
    intercept = 1 - (slope*(GL(1)));
    SI = round(imlincomb(slope,I,intercept,'double'));
end
SI(SI > NL) = NL;
SI(SI < 1) = 1;
numOffsets = size(Offset,1);

if NL ~= 0
    s = size(I);
    [r,c] = meshgrid(1:s(1),1:s(2));
    r = r(:);
    c = c(:);
    % Compute GLCMS
    GLCMS = zeros(size(Iorig,1),size(Iorig,2),NL,NL,numOffsets,'single');
    for k = 1 : numOffsets
%         fprintf('========%d=======',k);
        GLCMS(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),:,:,k) = computeGLCM(r,c,Offset(k,:),SI,NL,NeighborSize);
%         squeeze(GLCMS(InsideRange(1,1)+55,InsideRange(2,1)+15,:,:,k))
    end
else
	GLCMS = zeros(0,0,numOffsets);
end

end
%--------------------------------------------------------------------------
function [oneGLCM] = computeGLCM(r,c,offset,si,nl,NeighborSize)
oneGLCM = zeros(size(si,1),size(si,2),nl,nl);
Neigh = ones(NeighborSize*2+1);
r2 = r + offset(1);
c2 = c + offset(2);
[nRow nCol] = size(si);
outsideBounds = find(c2 < 1 | c2 > nCol | r2 < 1 | r2 > nRow);
Index1 = r + (c - 1)*nRow;
Index1(outsideBounds) = [];
r2(outsideBounds) = []; 
c2(outsideBounds) = [];
Index2 = r2 + (c2 - 1)*nRow;
hasNoValue = isnan(si(Index1)) | isnan(si(Index2));
Index1 = Index1(~hasNoValue);
Index2 = Index2(~hasNoValue);
for i=1:nl
    for j=1:nl
        hasVec = zeros(size(si));
        Index = Index1(si(Index1)==i&si(Index2)==j);
        if(isempty(Index))
            oneGLCM(:,:,i,j) = zeros(size(si));
        else
            hasVec(Index) = 1;
            oneGLCM(:,:,i,j) = conv2(hasVec,Neigh,'same');
        end
    end
end
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
