function qubToTxt( filename )
% Convert forQuB file to an auto.txt format file for use with makeplots,
% etc

if ~exist('filename','var')
    f = getFiles;
    filename = f{1};
end

SAVETRACES( FNAME, 'txt',    D,A,F,IDs )

fret = load(filename);

f = inputdlg('What is the trace length (in frames) for this data?');
traceLen = str2double(f)

nTraces = length(fret)/traceLen;

fret = reshape(f, nTraces, traceLen);

[p,name] = fileparts(filename);
idx = 1:nTraces;
x = cell(nTraces,1);
x = strcat(name,x);
ids = strcat(x,num2str(idx','%1.0f'));

d = zeros(size(fret));
a = zeros(size(fret));


outputFilename = strrep(filename, '.qub.txt','.txt');

saveTraces( outputFilename, d,a,fret,ids );

end
