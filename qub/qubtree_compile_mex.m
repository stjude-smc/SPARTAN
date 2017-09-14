% Compile qub_loadModel and qub_saveModel mex functions.
% These are statically linked to the QubTree library (see qubtree directory).
%
% To make the qubtree libraries, see qubtree/Makefile for Linux/MacOSX and the
% associated VisualStudio project files for Windows. The libraries will be
% created in the same directory as the source code (qubtree/).
%
% -largeArrayDims forces 64-bit addressing support, which is included for
% future compatibility (this setting will be default soon).

%   Copyright 2007-2017 Cornell University All Rights Reserved.



if strcmp(mexext,'mexw64')
    % Windows 64-bit
    mex -O -win64 -largeArrayDims -Lqubsuite\ -lqubtree -outdir ..\binary\ qub_loadTree.cpp treestruct.cpp
    mex -O -win64 -largeArrayDims -Lqubsuite\ -lqubtree -outdir ..\binary\ qub_saveTree.cpp treestruct.cpp

elseif strcmp(mexext,'mexa64') || strcmp(mexext,'mexmaci64')
    % 64-bit Linux
    mex -O GCC='/usr/bin/gcc-4.7' -largeArrayDims -Lqubsuite/ -lqubtree -outdir ../binary qub_loadTree.cpp treestruct.cpp
    mex -O GCC='/usr/bin/gcc-4.7' -largeArrayDims -Lqubsuite/ -lqubtree -outdir ../binary qub_saveTree.cpp treestruct.cpp

elseif strcmp(mexext,'mexa64') || strcmp(mexext,'mexmaci64')
    % 64-bit Mac (Intel)
    mex -O -largeArrayDims -Lqubsuite/ -lqubtree -outdir ../binary qub_loadTree.cpp treestruct.cpp
    mex -O -largeArrayDims -Lqubsuite/ -lqubtree -outdir ../binary qub_saveTree.cpp treestruct.cpp

else
    disp('Failed: unsupported architecture for compiling');
end





