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
    
end %public methods

end %class Traces






