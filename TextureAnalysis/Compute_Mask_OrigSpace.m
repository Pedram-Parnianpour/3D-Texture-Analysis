function Compute_Mask_OrigSpace(varargin)
invDeformMapFiles = ParseParams(varargin);
batchName = 'data/batch_Mask.mat';
Thresh = 0.1;
batchFile = batchName;
BATCH = load(batchFile);
spm_jobman('initcfg');
fprintf('Starting Time: %s\n',datestr(now));
for i=1:length(invDeformMapFiles)
    [outpath outName] = fileparts(invDeformMapFiles{i});
    invDeformMap = [outpath '/' outName '.nii'];
    outName = strrep(outName, 'iy_r', '');
    maskPathName = [outpath '/' outName '_Mask.nii'];
    BATCH.matlabbatch{1}.spm.util.defs.ofname = outName;
    BATCH.matlabbatch{1}.spm.util.defs.comp{1}.def{1} = invDeformMap;
    spm_jobman('run',BATCH.matlabbatch);
    currentDir = pwd;
    movefile([currentDir '/wT1_Mask.nii'],maskPathName);
    delete([currentDir '/y_' outName '.nii']);
    vol = load_untouch_nii(maskPathName);
    img = vol.img>Thresh;
    [x,y,z] = ndgrid(-3:3);
    se = strel(sqrt(x.^2 + y.^2 + z.^2) <=3);
    img = imclose(img,se);
    img = imopen(img,se);
    % img = imfill(img,'holes');
    vol.img = img;
    vol.hdr.dime.glmax = max(img(:));
    vol.hdr.dime.glmin = min(img(:));
    vol.hdr.hist.descrip = 'Mask';
    save_untouch_nii(vol,maskPathName);    
end
fprintf('End Time: %s\n',datestr(now));
end
%--------------------------------------------------------------------------
function invDeformMapFiles = ParseParams(varargin)
invDeformMapFiles = cell(1,length(varargin{1}{1}.invdefmap));
for i=1:length(invDeformMapFiles)
    invDeformMapFiles{i} = strtok(varargin{1}{1}.invdefmap{i},',');
end
end



