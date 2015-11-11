classdef TracesFret < Traces
% TracesFret: two-color FRET traces
% 
%    Traces class for two-color FRET experiments, with data channels for donor
%    and acceptor fluorescence and FRET efficiency (fret). These data channels
%    may be empty if any are not used.
%
%    A new traces object, with all data filled with zeros, can be created as
%    follows, where the default channels are donor, acceptor, and fret.
% 
%          traces = Traces(nTraces,nFrames);
%
%    If only a subset of fields is valid, add an argument for channel names:
%
%          traces = Traces(nTraces,nFrames,channelNames);
%
%    If corrections have been made to the fluorescence data and you want to
%    recalculate FRET values, use this method. Thresholds of total intensity
%    below which FRET is defined as zero can be given as a parameter (optional).
%
%         data.recalculateFret( thresholds );
%
%    See Traces.m for more information on Traces objects.

%   Copyright 2007-2015 Cornell University All Rights Reserved.


% Declaration and initialization of public parameters
properties (SetAccess=public, GetAccess=public)
    % Common properties that almost all data will have.
    donor    = [];
    acceptor = [];
    fret     = [];
end %end public properties


% Properties calculated as needed, but not actually saved in the class.
properties (Dependent)
    idxFluor;  %indexes to fluorescence (not FRET) channels)
    idxFret;   %indexes to FRET channels
end



methods
    %% ================ CONSTRUCTORS ================ %%
    function this = TracesFret(varargin)
        % Constructors cannot create an object of a different class than their
        % own. Instead the data are copied into an instance of this class.
        if nargin>=1 && (ischar(varargin{1}) || isa(varargin{1},'Traces')),
            args = {};  %will end up creating an empty instance
        else
            args = varargin;
            
            % Assign default parameter values
            if numel(args)<2,
                args = {0,0};
            end

            if numel(args)<3,
                args{3} = {'donor','acceptor','fret'};
            end
        end
        
        % Call superclass constructor to set up the object.
        this = this@Traces( args{:} );
        
        % -------- Load traces from file --------
        if nargin==1 && ischar(varargin{1}),
            this = loadTraces(varargin{1});
            
        % -------- Copy Constructors --------
        % This will raise an error if the user attempts to load a .traces file
        % with a different type than this class.
        elseif nargin==1 && strcmp(class(varargin{1}),class(this))
            this = copy(varargin{1});
            
        elseif nargin==1 && isa(varargin{1},'Traces')
            % FIXME: this is broken. But unused?
            this = copy(varargin{1});
        end
    end
    
    
    
     %% ================  GET/SET METHODS  ================ %%
     
    function idx = get.idxFluor(this)
        idx = find(  cellfun( @(x) isempty(strfind(x,'fret')), this.channelNames )  );
    end
    
    function idx = get.idxFret(this)
        idx = find(  cellfun( @(x) ~isempty(strfind(x,'fret')), this.channelNames )  );
    end

    % Calculate total fluorescence intensity
    function T = total(this,varargin)
        % Create struct for subsref, if the user requested a slice of the total
        % intensity matrix, saving memory by not constructing the full matrix.
        % The syntax is typically: T = data.total(1:10,:)
        S(1).type = '.';
        if nargin>1,
            S(2).type = '()';
            S(2).subs = varargin;
        end
        
        % Sum intensity from all fluorescence channels.
        for c=1:numel(this.idxFluor),
            S(1).subs = this.channelNames{this.idxFluor(c)};
            if c==1,
                T = subsref( this, S );
            else
                T = T + subsref( this, S );
            end
        end
    end
    
    
    %% ================ DATA MANIPULATION ================ %%
    
    % Add trace data to the current object, modifying it in place.
    function this = appendTraces( this, donor, acceptor, fret, traceMetadata )
        % Make sure all inputs are valid and compatible
        checkValid(this);
        assert( all([size(donor,2),size(acceptor,2),size(fret,2)]==this.nFrames) & ...
                all([size(acceptor,1),size(fret,1)]==size(donor,1)), 'Trace size mismatch' );
        
        % Append new data fields
        this.donor    = vertcat(this.donor,donor);
        this.acceptor = vertcat(this.acceptor,acceptor);
        this.fret     = vertcat(this.fret,fret);
        
        % Append new metadata fields, if possible.
        if nargin>=5,
            if numel(this.traceMetadata)~=size(donor,1),
                warning('Trace metadata size mismatch. Ignoring');
            elseif ~isempty(setdiff( fieldnames(this.traceMetadata), traceMetadata )),
                % FIXME: allow partial metadata list of metadata fields.
                warning('Trace metadata field mismatch. Ignoring' );
            else
                this.traceMetadata = vertcat( to_col(this.traceMetadata), to_col(traceMetadata) );
            end
        end
        
        % If traceMetadata is missing values, pad the end with empty values.
        nPad = this.nTraces - numel(this.traceMetadata);
        if nPad > 0,
            [this.traceMetadata(end+1:end+nPad).ids] = deal('');
            if ~iscolumn(this.traceMetadata),
                this.traceMetadata = reshape(this.traceMetadata,[this.nTraces 1]);
            end
        end
        
        checkValid(this);
    end
    
    
    function this = recalculateFret( this, idx, varargin )
    % DATA.recalculateFret() recalculates FRET efficiency from fluorescence
    % in all traces in the TracesFret object DATA. FRET to zero during donor
    % blinking events using automatically-calculated thresholds.
    % 
    % DATA.recalculateFret(IDX) only alters the subset of traces IDX.
    % 
    % DATA.recalculateFret(..., THRESH) recalculates FRET effciency using 
    % the supplied thresholds (used by sorttraces).
    % 
    % The method used for determining where the donor is dark is specified
    % in fileMetadata.zeroMethod, which can be 'skm' or 'threshold' (default).
    
        if ~isChannel(this,'fret'),
            return; %no FRET to recalculate (single-color)
        end
    
        % Get list of traces to correct; correct all if not specified.
        if nargin<3,  
            idx = true(this.nTraces,1);
        end
        if isempty(idx), return; end
        assert( isvector(idx) );
        
        % Realculate FRET efficiency, only in the selected traces.
        this.fret(idx,:) = this.acceptor(idx,:)./this.total(idx,:);
        this.fret( isnan(this.fret) ) = 0;  %NaN can happen in low SNR traces.
        
        % Set FRET to zero when the donor is dark (total intensity at baseline).        
        if isfield(this.fileMetadata,'zeroMethod') && strcmpi(this.fileMetadata.zeroMethod,'skm')
            alive = skmTotal( this.total(idx,:), varargin{:} );
        else
            alive = thresholdTotal( this.total(idx,:), varargin{:} );
        end
        
        mask = false( size(this.fret) );
        mask(idx,:) = ~alive;
        this.fret(mask) = 0;
        
    end %METHOD recalculateFret
    

end %methods

end %classdef


