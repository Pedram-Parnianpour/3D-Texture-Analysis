function GLCM_ROI(varargin)

[OutDirectory,fileList,maskFile,GLCMOPT] = ParseParams(varargin);
Onemask = false;
if(length(maskFile)==1)
    MaskVol = load_untouch_nii(maskFile{1});
    mask = MaskVol.img>0;
    [d1,d2,d3] = ind2sub(size(mask), find(mask));
    InsideRange = [min(d1) max(d1);min(d2) max(d2);min(d3) max(d3)];
    mask = mask(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));    
    Onemask = true;
end

fprintf('Starting Time: %s\n',datestr(now));
Result = zeros(length(fileList),length(GLCMOPT.glcm_properties),'single');
outName = cell(length(fileList),1);
for i=1:length(fileList)
    if(~Onemask)
        MaskVol = load_untouch_nii(maskFile{i});
        mask = MaskVol.img>0;        
        [d1,d2,d3] = ind2sub(size(mask), find(mask));
        InsideRange = [min(d1) max(d1);min(d2) max(d2);min(d3) max(d3)];
        mask = mask(InsideRange(1,1):InsideRange(1,2),InsideRange(2,1):InsideRange(2,2),InsideRange(3,1):InsideRange(3,2));
    end
    vol = load_untouch_nii(fileList{i});
    vol = single(vol.img);
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
    newLower = (1/(GLCMOPT.quantLevel-1));
    newUpper = 1-newLower;
    I = (I-norm_fact_L)/(norm_fact_U-norm_fact_L)*(newUpper-newLower)+newLower;
    I(I>1) = 1;
    I(I<0) = 0;
    vol(mask) = I;
    vol(~mask(:)) = 0;    
    GLCM = CreateGLCM(vol,GLCMOPT.quantLevel,[0 1],GLCMOPT.D,mask);    
    texture = computeFeature(GLCM,GLCMOPT.glcm_properties);  
    Result(i,:) = texture;
    time_total = toc(timerStart);
    fprintf('Total time: %2.2f\n',time_total);
    [~,outName{i}] = fileparts(fileList{i});  
end
fprintf('Writing the results...\n');
% Writing the results
xlsFile = [OutDirectory 'ROI_D' num2str(GLCMOPT.D) 'ROI_Q' num2str(GLCMOPT.quantLevel)];
xlsCols = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
xlswrite(xlsFile,GLCMOPT.glcm_properties,['B1:' xlsCols{length(GLCMOPT.glcm_properties)+1} num2str(1)]);    
xlswrite(xlsFile,outName,['A2:A' num2str(length(fileList)+1)]);    
xlswrite(xlsFile,Result,['B2:' xlsCols{length(GLCMOPT.glcm_properties)+1} num2str(length(fileList)+1)]);    
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
function glcm = CreateGLCM(I, NL, GL,D,mask)
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
clear I SI;
% tic
[OffsetsMatrixS OffsetsMatrixE] = AllOffsets(vol_Gray,D,mask);
% g=toc
InsideMask = OffsetsMatrixS>0&OffsetsMatrixE>0;
glcm = zeros(NL,NL,'single');
fprintf('Computing GLCM: ')
n=0;
for i=1:NL
    for j=1:NL
        G = zeros(size(InsideMask),'uint8');
        G(InsideMask) = vol_Gray(OffsetsMatrixS(InsideMask))==i&vol_Gray(OffsetsMatrixE(InsideMask))==j;
        glcm(i,j) = sum(G(:));
    end
    fprintf(repmat('\b',1,n));
    msg = sprintf('%2.2f',i/NL*100);
    fprintf([msg '%%']);
    n=numel(msg)+1;                    
end
fprintf('\n');
clear vol_Gray;
norm_fact = sum(glcm(:));
glcm = glcm./norm_fact;
end
%--------------------------------------------------------------------------
function [OffsetsMatrixS OffsetsMatrixE] = AllOffsets(I,D,mask)
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
% For symmetric glcm we need half of offsets
% Offsets = Offsets(1:size(Offsets,1)/2,:);
OffsetsMatrixS = single(repmat(IND_s,[1 ,size(Offsets,1)]));
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
    OffsetsMatrixE(isInMask,ofst) = IND_e(isInMask2);
    fprintf(repmat('\b',1,n));
    msg = sprintf('%2.2f',ofst/size(Offsets,1)*100);
    fprintf([msg '%%']);
    n=numel(msg)+1;                
end
fprintf('\n');
end
%--------------------------------------------------------------------------
function [outDir,fileList,maskFiles,GLCMOPT] = ParseParams(varargin)
outDir = varargin{1}{1}.dir{1};
fileList = cell(1,length(varargin{1}{1}.data));
for i=1:length(fileList)
    fileList{i} = strtok(varargin{1}{1}.data{i},',');
end
maskFiles = cell(1,length(varargin{1}{1}.masks));
for i=1:length(maskFiles)
    maskFiles{i} = strtok(varargin{1}{1}.masks{i},',');
end
if(length(maskFiles)~=1 && length(maskFiles)~=length(fileList))
    errordlg('The number of mask(s) should be 1 or the same as the number of volumes!');
    finish;
end
GLCMOPT.D = varargin{1}{1}.optsROI.distance;
GLCMOPT.quantLevel = varargin{1}{1}.optsROI.quantlevel;
glcm_properties = {'autoc','contr','corrm','corrp','cprom','cshad','dissi','energ','entro','homom','homop','maxpr','sosvh','savgh','svarh','senth','dvarh','denth','inf1h','inf2h','indnc','idmnc'};
GLCMOPT.glcm_properties = {};
r=1;
for i=1:length(glcm_properties)
    [~,present] = evalc(['varargin{1}{1}.features.' glcm_properties{i}]);
    if(present==1)
        GLCMOPT.glcm_properties{r} = glcm_properties{i};
        r = r+1;
    end
end

end
