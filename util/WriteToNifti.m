function [ nii ] = WriteToNifti( hdr, fpath, x, varargin )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    p = inputParser();
    addRequired(p, 'hdr');
    addRequired(p, 'fpath');
    addRequired(p, 'x');
    addOptional(p, 'ind', [], @isnumeric);
    addOptional(p, 'datatype', 'float', @ischar);
    parse(p, hdr, fpath, x, varargin{:});

    % Scaling is used to help express a range of numbers without truncating
    % when a data type does not have sufficient range or precision to capture
    % the raw data. For example, when representing naturally floating-point
    % data in an integer format. For now, I assume data will be written without
    % scaling.
    hdr.dime.scl_slope = 1;
    hdr.dime.scl_inter = 0;
    switch p.Results.datatype
        case 'binary'
            % 1 Binary (ubit1, bitpix=1)
            % DT_BINARY
            hdr.dime.datatype = 1;
            hdr.dime.bitpix = 1;
            x = logical(x);

        case 'uint8'
            % Unsigned char (uint8, bitpix=8; integers 0--255)
            % DT_UNSIGNED_CHAR, DT_UINT8
            hdr.dime.datatype = 2;
            hdr.dime.bitpix = 8;
            x = uint8(x);
            
        case 'int'
            % Signed integer (int32, bitpix=32; integers -2,147,483,648 -- 2,147,483,647)
            % DT_INT32, NIFTI_TYPE_INT32
            hdr.dime.datatype = 8;
            hdr.dime.bitpix = 32;
            x = int32(x);

        case 'float'
            % Floating point (single or float32, bitpix=32)
            % DT_FLOAT32, NIFTI_TYPE_FLOAT32
            hdr.dime.datatype = 16;
            hdr.dime.bitpix = 32;
            x = single(x);

        case 'double'
            % Floating point    (single or float32, bitpix=32)
            % DT_FLOAT32, NIFTI_TYPE_FLOAT32
            hdr.dime.datatype = 64;
            hdr.dime.bitpix = 64;
            x = double(x);

    end
    if nargin < 4 || isempty(p.Results.ind)
        X = x;
    else
        X = VectorTo3D(hdr, x, p.Results.ind);
    end
    nii = make_nii(X);
    nii.hdr = hdr;
    save_nii(nii,fpath);
end

