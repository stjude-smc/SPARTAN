classdef MovieParser < handle
% MovieParser   Display an interactive list of fret traces.
%
%   p = MovieParser(FILE) opens the given fluorescence movie file (.tif or .stk)
%   and prepares a parser object that will assist in movie analsis tasks, such
%   as summing fluorescence fields, detecting molecule locations, finding
%   integration windows for each, and saving the traces to file.
%
%   See also: Movie, ChannelExtractor, gettraces_gui.

%   Copyright 2017-2022 All Rights Reserved.



properties (GetAccess=public, SetAccess=protected)
    % Basic data available when movie is first loaded in openStk().
    chExtractor;         % ChannelExtractor encapsulates Movie and splits frame data into channels
    stdbg;               % stdev of background noise at the end of movie
    
    % Picked molecules from getPeaks()
    total_t;             % Registered, total intensity image used for picking
    peaks;               % Molecule locations in all fields (molID,dim,channel)
    total_peaks;         % ... in total intensity image
    rejectedPicks;       % Locations of molecules with overlapping PSFs
    rejectedTotalPicks;  % ... in total intensity image
    fractionOverlapped;  % Fraction of molecules rejected due to overlapping PSFs
    alignStatus;         % Alignment struct: dx, dy, theta, sx, sy, abs_dev, tform, quality
    
    % Integration windows from getIntegrationWindows()
    regionIdx;           % X,Y coordinates of each integration window
    bgMask;              % Logical mask of pixels used for summing background fluorescence
    
    % Integration startistics for display, also from getIntegrationWindows()
    integrationEfficiency;  % Estimated fraction of intensity collected
    psfWidth;               % Average number of pixels to integrate 70% of total intensity
    fractionWinOverlap;     % Fraction of pixels used by multiple molecules
end


% These can be directly manipulated by the user (gettraces_gui.m)
% FIXME: use set method to update structure and raise event if changed.
properties (GetAccess=public, SetAccess=public)
    % General analysis settings. See cascadeConstants.m.
    % The values in params include all known cameras for a microscope,
    % not just the ones in the current file.
    params;
    
    % Each channel's assignment for trace interpretation
    roles;
    
    % Index into params of each channel actually in the current movie.
    idxActiveChannels;
end


properties (Dependent)
    % Number of frames in the movie after deinterlacing colors.
    % This will be equal to movie.nFrames unless fluorescence channels
    % appear as interlaced frames instead of field areas (not common).
    nFrames;
    
    nChannels;
end


properties (Constant)
    validRoles = {'donor','acceptor','donor2','acceptor2','factor','factor2'};
end



methods
    % Constructor
    function this = MovieParser(varargin)
        if nargin>0
            this = this.openStk( varargin{:} );
        end
    end
    
    % Open movie and prepare for viewing
    this = openStk(this, input, params);
    
    % Detect peaks of intensity in registered, total intensity image
    this = getPeaks(this);
    
    % Find integration windows for each molecule
    this = getIntegrationWindows(this);
    
    % Sum fluorescence in integration windows and save fluorescence traces
    integrateAndSave(this, filename);
    
    % get/set methods
    function value = get.nFrames(this)
        % Accounts for channels stacked as separate frames
        value = this.chExtractor.nFrames;
    end
    
    function value = get.nChannels(this)
        % Accounts for channels stacked as separate frames
        value = numel(this.idxActiveChannels);
    end
    
    
    
    %% Update correction parameters
    
    function updateCrosstalk(this)
    % Prompt user for spectral crosstalk values.

        idxCh = this.idxActiveChannels;
        if this.nChannels<2, return; end  %crosstalk not defined for single color

        % Extract crosstalk sub-matrix associated with active channels
        crosstalk = this.params.crosstalk( idxCh, idxCh );

        % Enumerate all possible (forward) crosstalk pairs.
        [src,dst] = find(  triu(true(this.nChannels),1)  );

        % Prompt the user for crosstalk parameters.
        % Channel names MUST be in order of wavelength!
        prompts  = cell(numel(src), 1);
        defaults = cell(numel(src), 1);
        names = {this.chExtractor.channels.name};

        for i=1:numel(src),
            prompts{i}  = sprintf( '%s to %s', names{src(i)}, names{dst(i)} );
            defaults{i} = num2str( crosstalk(src(i),dst(i)) );
        end

        result = inputdlg(prompts, 'gettraces', 1, defaults);
        if isempty(result), return; end  %user hit cancel
        result = cellfun(@str2double, result);

        % Verify the inputs are valid and save them.
        if any( isnan(result) | result>1 | result<0 ),
            warndlg('Invalid crosstalk values');
            return;
        end

        for i=1:numel(src),
            crosstalk(src(i), dst(i)) = result(i);
        end
        this.params.crosstalk(idxCh,idxCh) = crosstalk;

    end  %function updateCrosstalk
    
    
    function updateScaling(this)
        % Set values for scaling the fluorescence intensity of acceptor channels so
        % that they all have the same apparent brightness (gamma=1).

        idxCh = this.idxActiveChannels;

        % Extract subset of channels in the current movie.
        scaleFluor = this.params.scaleFluor(idxCh);

        % Prompt the user for new multipliers for gamma correction.
        defaults = cellfun(@num2str, num2cell(scaleFluor), 'UniformOutput',false);
        result = inputdlg( {this.chExtractor.channels.name}, 'gettraces', 1, defaults );
        if isempty(result), return; end

        % Verify the inputs are valid and save them.
        result = cellfun(@str2double, result);
        if any(isnan(result)),
            warndlg('Invalid scaling values');
            return;
        end

        this.params.scaleFluor(idxCh) = to_row(result);

    end %function updateScaling

    
    
    %% Update ChannelExtractor parameters
    
    % Prompt user to alter parameter values for a specific channel
    function updateChannel(this, chID)
        % Collect current channel parameters as prompt defaults
        idxParams = this.idxActiveChannels(chID);
        channel = this.params.channels(idxParams);

        channel.role = this.roles{chID};
        fields = {'role','name','wavelength','photonsPerCount'};
        prompt = {'Role:', 'Name:', 'Wavelength (nm):', 'photonsPerCount'};
        types  = { this.validRoles, {this.params.channels.name}, @isnumeric, @isnumeric };
        channel = settingdlg(channel, fields, prompt, types, 'Field Settings');
        
        % Save the new parameter values
        if ~isempty(channel)
            this.roles{chID} = channel.role;
            channel = rmfield(channel,'role');
            this.params.channels(idxParams) = channel;
            this.chExtractor.channels(chID) = channel;
        end
    end
    
    % Prompt user for field arrangement and channel assignments
    function success = updateFieldArrangement(this)
        
        % Replace field arrangement matrix with indexes into the full list
        % of allowed channels.
        geo = this.chExtractor.fieldArrangement;
        profileNames = {this.params.channels.name};
        
        for i=1:numel(geo)
            if geo(i)==0, continue; end
            chName = this.chExtractor.channels(geo(i)).name;
            geo(i) = find(  strcmp(chName, profileNames)  );
        end
        
        % Prompt user for field geometry and channel names.
        geo = fieldArrangementDialog( geo, {this.params.channels.name} );
        success = ~isempty(geo);
        if ~success, return; end
        
        % Reset channelExtractor and roles.
        this.idxActiveChannels = sort( geo(geo>0) );
        this.chExtractor.channels = this.params.channels( this.idxActiveChannels );
        
        % Replace field arrangement matrix with indexes into only the
        % channels actually used.
        geoUsed = sort( geo(geo>0) );
        for i=1:numel(geo)
            if geo(i)==0, continue; end
            geo(i) = find( geo(i)==geoUsed );
        end
        this.chExtractor.fieldArrangement = geo;
        
        % Make a guess as to the new roles
        nCh = numel( this.idxActiveChannels );
        if nCh<4
            temp = {'donor','acceptor','acceptor2'};
            this.roles = temp(1:nCh);
        else
            this.roles = {'donor','acceptor','donor2','acceptor2'};
        end
        
        %disp(geo);  disp( {this.chExtractor.channels.name} );

        this.chExtractor.avgTop( this.params.nAvgFrames, this.params.subtractBGImage );
    end


    function updateAlex(this)
    % Open dialog to alter illumination metadata (for ALEX experiments)
        a = inputdlg('Enter wavelength series in one cycle:','gettraces');
        if ~isempty(a)
            a = cellfun( @str2double, strsplit(a{1},{',',' '}) );
            if any(isnan(a)) && numel(a)==numel(unique(a))
                errordlg('Invalid input. Must be comma-separated list of unique wavelengths');
                return;
            end
            
            % Use wavelength list to build illumination metadata fields
            lasers = struct('wavelength',num2cell(sort(a)));
            for i=1:numel(lasers)
                firstFrame = find( lasers(i).wavelength==a );
                lasers(i).framesActive = firstFrame:numel(lasers):this.chExtractor.nFrames;
            end
            this.chExtractor.lasers = lasers;
        end
    end
end



end %class MovieParser


