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
% 



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
            this.copyDataFrom(varargin{1});
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
    function T = total(this)
        T = zeros(this.nTraces,this.nFrames);
        for c=1:numel(this.idxFluor),
            T = T + this.( this.channelNames{this.idxFluor(c)} );
        end
    end
    
    
    %% ================ DATA MANIPULATION ================ %%
    % These are shortcut functions for common tasks of manipulating fluorescence
    % and FRET data, in particular to make corrections for crosstalk and gamma.
    % The data stored in the Traces object is always corrected, and we assume
    % that any correcctions are stored in traceMetadata correctly, but this
    % might not be the case... FIXME: is there anything we can do about that??
    
    function this = recalculateFret( this, thresholds )
    % Recalculate FRET efficiencies (all fields) using fluorescence. 
    %    data.recalculateFret( THRESH );
    % Replaces the current "fret" property with the newly calculated traces.
    % Thresholds of total intensity below which FRET is set to zero can be given
    % as the extra parameter (THRESH, Nx1 array, one per trace). Otherwise, they
    % calculated from the data as a few standard deviations above background.
    %
    % FIXME: would this require that fretThreshold (manually tuned) be included
    % in traceMetadata??? Is that value even relevant after doing the
    % corrections?
    
        if ~isChannel(this,'fret'), warning('Not FRET data?'); end
        assert( ~isChannel(this,'acceptor2') & ~isChannel(this,'donor2'), ...
                'Not valid for 3-color data (yet, FIXME)' );
        
        constants = cascadeConstants;
        
        if nargin<2,
            thresholds = zeros(data.nTraces,1);
        end
        
        % Determine the end of each trace (donor bleaching event.
        total = this.donor+this.acceptor;
        lt = max(1, calcLifetime(total) );
        
        % Calculate FRET for each trace.
        this.fret = this.acceptor./total;
        
        for i=1:this.nTraces,
            % Set FRET to zero after donor photobleaching.
            this.fret(i, lt:end ) = 0;
            
            % Set FRET to zero in areas where the donor is dark (blinking).
            % ie, when FRET is below a calculated threshold = 4*std(background)
            % This was taken from correctTraces.m
            % FIXME: what happens to manual FRET threshold adjustments??
            %        Should they be saved in traceMetadata?
            s = lt(i)+5;
            range = s:min(s+constants.NBK,this.nFrames);
            if numel(range)>=10,
                if nargin<2, 
                    thresholds(i) = constants.blink_nstd*std(total(i,range));
                end
                darkRange = total( i, 1:lt(i) ) <= thresholds(i);
                this.fret(i,darkRange) = 0;
            end            
        end
    end
    
    
    
end %public methods

end %class Traces






