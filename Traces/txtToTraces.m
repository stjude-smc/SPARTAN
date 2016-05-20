function txtToTraces( inputFilenames, outputFilenames )
% Converts old format traces to the new binary format.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Get input and output filenames
if nargin<1,
    inputFilenames = getFiles('*.txt','Select .txt traces files to convert.');
end
if numel(inputFilenames)==0, return; end

if nargin<2,
    outputFilenames = strrep(inputFilenames,'.txt','.traces');
end

assert( numel(inputFilenames)==numel(inputFilenames), 'Invalid number of output filenames' );


% Load data into memory and resave in new format.
for i=1:numel(inputFilenames),
    saveTraces( outputFilenames{i}, loadTraces(inputFilenames{i}) );
end

end
