function forQuB( files, fretCh )
%FORQUB   Saves FRET data in the .qub.txt format for importing into QuB
%
%   forQuB( FILES ) loads each .traces file in the cell array FILES and saves
%   a corresponding .qub.txt file containing just the fret data that can be
%   loaded into QuB.
%
%   forQuB( FILES, CH ) saves the data channel CH (typically fret or fret2).

%   Copyright 2007-2015 Cornell University All Rights Reserved.


if nargin<2,
    fretCh = '';
end

% Get file names from user
if nargin<1,
    files = getFiles();
end
if ischar(files),
    files={files};
end
nFiles = numel(files);


for i=1:nFiles,
    data = loadTraces( files{i} );
    [p,f] = fileparts(files{i});
    
    % For multi-color imaging, ask the user which channel to save.
    % Only ask this once to simplify saving many files.
    if isempty(fretCh) && isChannel(data,'fret2'),
        fretNames = data.channelNames(data.idxFret);
        a = questdlg('Which FRET channel should be saved?', ...
                     'forQuB',fretNames{:},'Cancel', 'fret');
        if strcmp(a,'Cancel'), return; end
        fretCh = a;
    end
    
    if isChannel(data,'fret2'),
        f = [f '_' fretCh];  %#ok<AGROW>
        fret = data.(fretCh)';
    else
        fret = data.fret';
    end
    
    % Remove extreme values that can confuse QuB.
    fret( fret<-0.5 ) = -0.5;
    fret( fret>1 ) = 1;

    % Get an output filename. Channel name is appended for multi-color FRET.
    outfile = fullfile(p, [f '.qub.txt']);
    
    if nFiles==1,
        [f,p] = uiputfile(outfile, 'Save .qub.txt file');
        if ~ischar(f), return; end
        outfile = fullfile(p,f);
    end
    
    % Save the data to file
    fid = fopen(outfile, 'w');
    fprintf(fid, '%f\n', fret(:));
    fclose(fid);
    
    fprintf('Saved %d traces in %s\n', size(fret,2), outfile);
end

