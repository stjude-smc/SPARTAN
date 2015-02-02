classdef TracesFret4 < TracesFret
% TracesFret: multi-channel FRET traces
% 
%    Traces class with extra properties for additional data channels, including
%    donor2, acceptor2, and fret2 for three-color FRET and two-pair FRET
%    experiments, and factor/factor2 channels for miscellaneous fluorescence
%    channels (for example for dye-labeled factor binding to the ribosome).
%
%    These data channels may be empty if any are not used. When querying the
%    trace data, you can iterate over all channels using 'channelNames' or over
%    a specific type of data using idxFret and idxFluor.
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
%    See Traces.m and TracesFret.m for more information on Traces objects.
% 



% Declaration and initialization of public parameters
properties (Access=public)
    % Other fields for multi-color fret. These may be empty if they are not
    % used. To know what fields are "alive", check the channelNames list or find
    % which of these fields are not empty.
    donor2    = [];
    acceptor2 = [];
    fret2     = [];
    factor    = [];
    factor2   = [];
end %end public properties


properties (Dependent)
    idxFactor;  %indexes (in channelNames) of factor-binding channels.
end



methods
    %% ================ CONSTRUCTORS ================ %%
    function this = TracesFret4(varargin)
        % Constructors cannot create an object of a different class than their
        % own. Instead the data are copied into an instance of this class.
        if nargin>=1 && (ischar(varargin{1}) || isa(varargin{1},'Traces')),
            args = {0,0,{}};
        else
            args = varargin;
        end
        
        % Call superclass constructor to set up the object.
        this = this@TracesFret( args{:} );
                
        % -------- Load traces from file --------
        if nargin==1 && ischar(varargin{1}),
            this.copyDataFrom( loadTraces(varargin{1}) );
            
        % -------- Copy Constructors --------
        % This does not allow initializing TracesFret4 object with TracesFret
        % data, which seems like a valid thing to do! Doing so may require
        % copying property fields.
        elseif nargin==1 && strcmp(class(varargin{1}),class(this))
            this = copy(varargin{1});
            
        elseif nargin==1 && isa(varargin{1},'Traces')
            this.copyDataFrom(varargin{1});
        end
    end
    
    
     %% ================  GET/SET METHODS  ================ %%
     
     function idx = get.idxFactor(this)
        idx = find(  cellfun( @(x) ~isempty(strfind(x,'factor')), this.channelNames )  );
     end
     
    
    %% ================ DATA MANIPULATION ================ %%
    
    function this = recalculateFret( this, thresholds )
    % Recalculate FRET efficiencies (all fields) using fluorescence. 
    %    data.recalculateFret( THRESH );
    % Replaces the current "fret" properties with the newly calculated traces.
    % Thresholds of total intensity below which FRET is set to zero can be given
    % as the extra parameter (THRESH, Nx1 array, one per trace). Otherwise, they
    % calculated from the data as a few standard deviations above background.
    % FIXME: specify that THRESH is NxM, where M is the number of FRET channels.
    %
    % NOTE: for three-color FRET, there are two ways to calculate it:
    %   1) We assume there is no donor->acceptor2 crosstalk. Here we can
    %   actually calculate FRET efficiencies.
    %   2) Otherwise, we can't actually calculate FRET, but as a placeholder, we
    %   can calculate the fraction of total intensity from each acceptor.
    % This function will look in fileMetadata to find out which to use, and if
    % not known, will ask the user to clarify.
    % 
            
        % If there is only one FRET channel, just use the subclass method.
        if ~isChannel(this,'acceptor2'),
            recalculateFret@TracesFret( this, thresholds );
            return;
        end
        
        % If this is four-color FRET (two totally independent FRET pairs),
        % calculate FRET for each pair separately and then combine.
        % NOTE: this assumes independence! Use ALEX for if you can't assume it.
        if isChannel(this,'donor2'),
            error('STUB: recalculateFret with four-color FRET');
        end
        
        constants = cascadeConstants;
        
        if nargin<2,
            thresholds = zeros(this.nTraces,1);  %fret1 only
        else
            if size(thresholds,2)>1,
                warning('Only fret1 thresholds will be used');
            end
        end
        
        
        % Calculate FRET for each trace.
        assert( this.isChannel('fret') && this.isChannel('fret2') );
        acc = this.acceptor + this.acceptor2; %acceptors total intensity
        total = this.donor + acc;
        
        if isfield(this.fileMetadata,'isTandem3') && this.fileMetadata.isTandem3,
            % Actually FRET. Assumes no donor->acceptor2 FRET.
            this.fret  = acc./total;
            this.fret2 = this.acceptor2./acc;
            isTandem3 = true;
        else
            % Not FRET, instead calculate fraction of total intensity in each channel.
            this.fret  = this.acceptor./total;
            this.fret2 = this.acceptor2./total;
            isTandem3 = false;
        end
        
        
        % Set FRET to zero when donor is dark (below total intensity threshold).
        lt = max(1, calcLifetime(total) );  %point where donor bleaches in each trace.
        
        for i=1:this.nTraces,
            % Set FRET to zero after donor photobleachinsg.
            this.fret(i, lt(i):end ) = 0;
            this.fret2(i, lt(i):end ) = 0;
            
            % Set FRET to zero in areas where the donor is dark (blinking).
            % ie, when FRET is below a calculated threshold = 4*std(background)
            % FIXME: the threshold (if given) should be applied regardless.
            s = lt(i)+5;
            range = s:min(s+constants.NBK,this.nFrames);
            if numel(range)>=10 && nargin<2,
                thresholds(i) = constants.blink_nstd*std(total(i,range));
            end
            
            if thresholds(i)>0,  %if it could be found or was user-defined,                
                darkRange = total( i, 1:lt(i) ) <= thresholds(i);
                this.fret(i,darkRange) = 0;
                this.fret2(i,darkRange) = 0;
            end
            
            % FRET2 is zero when the the total acceptor signal is so low that
            % FRET can't be calculated correctly (such as after photobleaching).
            % This may work poorly if crosstalk correction is bad.
            % FIXME: these parameters should be defined in cascadeConstants.
            if isTandem3,
                this.fret2( i, this.fret(i,:)<0.2 ) = 0;

                fretRange = rleFilter( this.fret(i,:)>=0.25, constants.rle_min );
                fret2_end = find(fretRange,1,'last');
                if ~isempty(fret2_end)
                    this.fret2(i, fret2_end:end ) = 0;
                end
            end
        end
        
        % Remove any NaN values. These usually happen when calcLifetime can't
        % see the bleaching step (b/c of low SNR), calculating FRET values from
        % zero-intensity noise. Can this be fixed?
        this.fret( isnan(this.fret) ) = 0;
        this.fret2( isnan(this.fret2) ) = 0;
        
        % TODO: save the thresholds that were used in traceMetadata.
        % FIXME: This would require allowing matrices in traceMetadata!
    end
end %public methods

end %class Traces






