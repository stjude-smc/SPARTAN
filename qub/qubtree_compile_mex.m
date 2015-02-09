% Script to compile qub_loadModel and qubSave model mex functions. These
% are linked to a pre-compiled library of functions that are derived from
% the QuB source code. Where possible, use static libraries to minimize
% possible link failures.
%
% To make the qubtree libraries, see qubtree/Makefile for Linux/MacOSX and the
% associated VisualStudio project files for Windows. The libraries will be
% created in the same directory as the source code (qubtree/).
%
% -largeArrayDims forces 64-bit addressing support, which is included for
% future compatibility (this setting will be default soon).
%

if strcmp(mexext,'mexw64'),
    % Windows 64-bit, DLL
    mex -v -O -win64 -Lqubtree\ -lqubtree_w64 -outdir ..\binary qub_loadTree.cpp treestruct.cpp
    mex -v -O -win64 -Lqubtree\ -lqubtree_w64 -outdir ..\binary qub_saveTree.cpp treestruct.cpp


elseif strcmp(mexext,'mexw64'),
    % Mac OSX, Intel 64-bit, static
    mex -v -O -largeArrayDims -Lqubtree/ -lqubtree_static -outdir ../binary qub_loadTree.cpp treestruct.cpp
    mex -v -O -largeArrayDims -Lqubtree/ -lqubtree_static -outdir ../binary qub_saveTree.cpp treestruct.cpp
    
elseif strcmp(mexext,'mexa64'),
    % Linux, Intel 64-bit, static
    mex -v -O -largeArrayDims -Lqubtree/ -lqubtree_static -outdir ../binary qub_loadTree.cpp treestruct.cpp
    mex -v -O -largeArrayDims -Lqubtree/ -lqubtree_static -outdir ../binary qub_saveTree.cpp treestruct.cpp
    
else
    disp('Failed: unsupported architecture for compiling');
end





