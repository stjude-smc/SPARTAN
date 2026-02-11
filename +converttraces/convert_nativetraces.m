% conversion of .traces to .mat
[f,p] = uigetfile('*.traces', 'Select a traces file to analyze');
if isequal(f,0)
    error('No file selected. Cannot continue.');
end

inputFile = fullfile(p,f);
tracesToMat( inputFile)
