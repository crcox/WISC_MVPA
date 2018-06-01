function HTCondor_struct2nii(ResultStruct, NiftiHeaders, VariableToMap, varargin)
% HTCondor_struct2nii 
    p = inputParser();
    addRequired(p,'ResultStruct', @isstruct);
    addRequired(p,'NiftiHeaders', @isstruct);
    addRequired(p,'VariableToMap',@ischar);
    addParameter(p,'outdir','SolutionMaps',@ischar);
    addParameter(p,'filestring','%02d_%02d.nii',@ischar);
    addParameter(p,'filevars',{'subject','cvholdout'},@iscellstr);
    parse(p,ResultStruct, NiftiHeaders, VariableToMap, varargin{:});
    
    varToMap = p.Results.VariableToMap;
    if ~(exist(p.Results.outdir,'dir')==7)
        mkdir(p.Results.outdir);
    end
    
    for i = 1:numel(ResultStruct)
        R = ResultStruct(i);
        x = cell(numel(p.Results.filevars),1);
        for j = 1:numel(x)
            x{j} = R.(p.Results.filevars{j});
        end
        fname = sprintf(p.Results.filestring, x{:});
%         if isstruct(R.subject)
%             subject = R.subject.subject;
%             cvholdout = R.subject.cvholdout;
%         else
%             subject = R.subject;
%             cvholdout = R.cvholdout;
%         end
%         fname = sprintf(p.Results.filestring, ...
%             subject, ...
%             cvholdout);
        fpath = fullfile(p.Results.outdir,fname);
        hdr = NiftiHeaders(R.subject).hdr;
        hdr.dime.scl_slope = 1;
        hdr.dime.datatype = 16;
        hdr.dime.bitpix = 32;
        X = zeros(hdr.dime.dim(2:4));
        X(R.coords.ind) = R.(varToMap);
        nii = make_nii(X);
        nii.hdr = hdr;
        save_nii(nii,fpath);
    end
end 
