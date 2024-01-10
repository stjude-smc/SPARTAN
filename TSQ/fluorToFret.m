function fluorToFret( files, chName, normalization )
% cy5forQuB    Normalize fluorescence data and save in fret field.
%
%   Fills fret data field with intensity values that are normalized to the
%   mean total fluorescence intensity all traces before bleaching.
%
%   Copyright 2007-2022 All Rights Reserved.


%% Process input arguments

% Get filename from user if not specified
if nargin < 1,
    files = getFiles;
end
if ~iscell(files),
    files = {files};
end
nFiles = numel(files);


%% Perform data manipulation

for i=1:nFiles
    
    % Load FRET data file
    filename = files{i};
    data = loadTraces( filename );
    
    % If not obvious, prompt user which fluorescence field should be used.
    if i==1 && nargin<2
        if data.nChannels==1
            chName = 'total';
        else
            fname = [data.channelNames 'total'];
            index = listdlg('ListString',fname, 'SelectionMode','single');
            if isempty(index), return; end  %user hit cancel.
            chName = fname{index};
        end
    end
    
    % Normalize the intensity to the ensemble average.
    if nargin<3
        % Average of the first 10 frames, one value per trace.
        val = data.(chName)(:,1:10);
        normalization = mean( val, 2 );
        
        % Average of the first 10 frames, one value for the whole file.
        %val = data.(chName)(:,1:10);
        %normalization = mean( val(:) );
        
        % Mean total intensity (like autotrace): one value per file.
        %stats = traceStat( data );
        %normalization = mean([stats.t]);
        
        % Mean total intensity (like autotrace): one value per trace.
        %stats = traceStat( data );
        %normalization = [stats.t];
    end
    
    if ~data.isChannel('fret')
        data.addChannel( 'fret', data.(chName) ./ normalization );
    else
        data.fret = data.(chName) ./ normalization;
    end
    %disp(normalization);

    % Save data to file.
    [p,f,e] = fileparts(filename);
    outname = fullfile(p,[f '_' mfilename e]);
    saveTraces(outname,data);
end

% END FUNCTION

