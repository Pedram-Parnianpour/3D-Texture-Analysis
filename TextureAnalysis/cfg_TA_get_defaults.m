function varargout = cfg_TA_get_defaults(defstr)
% Get dafual values of TA
currentPath = fileparts(which(mfilename));
persistent defaults;
if isempty(defaults)
    defaults.space.def = 0;
    defaults.quantlevel.def = 8;
    defaults.neighborsize.def = 1;
    defaults.distance.def = 1;
    defaults.smoothing.def = 0;
    defaults.autoc.def = 1;
    defaults.contr.def = 1;
    defaults.dissi.def = 1;
    defaults.energ.def = 1;
    defaults.entro.def = 1;
    defaults.homom.def = 1;
    defaults.homop.def = 1;
    defaults.sosvh.def = 1;
    defaults.indnc.def = 1;
    defaults.idmnc.def = 1;
    defaults.maxpr.def = 1;
    defaults.savgh.def = 1;
    defaults.senth.def = 1;
    defaults.svarh.def = 1;
    defaults.denth.def = 1;
    defaults.dvarh.def = 1;
    defaults.corrm.def = 1;
    defaults.corrp.def = 1;
    defaults.cprom.def = 1;
    defaults.cshad.def = 1;
    defaults.inf1h.def = 1;
    defaults.inf2h.def = 1;
    defaults.mask.def = {[currentPath 'data/T1_Mask.nii']}; 
end
[~,varargout{1}] = evalc(['defaults.' defstr ';']);
end

