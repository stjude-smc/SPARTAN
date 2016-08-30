function output = meanFret(files)

if nargin<1, files = getFiles; end
if ischar(files), files={files}; end
if numel(files)<1, return; end
nFiles = numel(files);

output = zeros(nFiles,1);

for i=1:numel(files)
    data = loadTraces(files{i});
    fret = data.fret(:,1:100);
    fret = fret(fret>0.2);
    output(i) = median(fret);
end

if nargout<1,
    figure;
    plot(output);
end

end