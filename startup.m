addpath core_dlag
addpath core_pcca
addpath core_fa
addpath plotting
addpath post-selection_inference
addpath simulation
addpath simulation/dlag
addpath simulation/pcca
addpath util
addpath util/precomp
addpath util/fileRetrieval
addpath util/fileRetrieval/dlag
addpath util/fileRetrieval/pcca
addpath util/fileRetrieval/fa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following code checks for the relevant MEX files (such as .mexa64
% or .mexglx, depending on the machine architecture), and it creates the
% mex file if it can not find the right one.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Toeplitz Inversion
path(path,'util/invToeplitz');
% Create the mex file if necessary.
if ~exist(sprintf('util/invToeplitz/invToeplitzFastZohar.%s',mexext),'file')
    try
        eval(sprintf('mex -outdir util/invToeplitz util/invToeplitz/invToeplitzFastZohar.c'));
        fprintf('NOTE: the relevant invToeplitz mex files were not found.  They have been created.\n');
    catch
        fprintf('NOTE: the relevant invToeplitz mex files were not found, and your machine failed to create them.\n');
        fprintf('      This usually means that you do not have the proper C/MEX compiler setup.\n');
        fprintf('      The code will still run identically, albeit slower (perhaps considerably).\n');
        fprintf('      Please read the README file, section Notes on the Use of C/MEX.\n');
    end
end

% Posterior Covariance Precomputation
path(path,'util/precomp/makePautoSumFast');
% Create the mex file if necessary.
if ~exist(sprintf('util/precomp/makePautoSumFast/makePautoSumFast.%s',mexext),'file')
    try
        eval(sprintf('mex -outdir util/precomp/makePautoSumFast/ util/precomp/makePautoSumFast/makePautoSumFast.c'));
        fprintf('NOTE: the relevant precomp mex files were not found.  They have been created.\n');
    catch
        fprintf('NOTE: the relevant precomp mex files were not found, and your machine failed to create them.\n');
        fprintf('      This usually means that you do not have the proper C/MEX compiler setup.\n');
        fprintf('      The code will still run identically, albeit slower (perhaps considerably).\n');
        fprintf('      Please read the README file, section Notes on the Use of C/MEX.\n');
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path(path,'util/invChol');
% Create the mex file if necessary.
if ~exist(sprintf('util/invChol/invChol_mex.%s',mexext),'file')
    try
        if strcmpi(computer, 'PCWIN64') || strcmpi(computer, 'GLNXA64') || ...
                strcmpi(computer, 'MACI64')
            disp('64-bit detected, compiling with -largeArrayDims flag...');
            mex ('-outdir', 'util/invChol',...
                'util/invChol/invChol_mex.c', ...                
                ['-I' matlabroot '/extern', '/examples', '/refbook'],...
                [matlabroot '/extern', '/examples', '/refbook', '/fort.c'], ...
                '-lmwlapack', '-largeArrayDims');
            
        else
            disp('32-bit detected, compiling without -largeArrayDims flag...');
            mex('-outdir', 'util/invChol',...
                'util/invChol/invChol_mex.c',...
                ['-I' matlabroot '/extern', '/examples', '/refbook'],...
                [matlabroot '/extern', '/examples', '/refbook', '/fort.c'], ...
                '-lmwlapack');
        end
        fprintf('NOTE: the relevant invChol mex files were not found.  They have been created.\n');
    catch
        fprintf('NOTE: the relevant invChol mex files were not found, and your machine failed to create them.\n');
        fprintf('      This usually means that you do not have the proper C/MEX compiler setup.\n');
        fprintf('      The code will still run identically, albeit slower (perhaps considerably).\n');
        fprintf('      Please read the README file, section Notes on the Use of C/MEX.\n');
    end
end