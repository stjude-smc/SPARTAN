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


properties (Constant,Hidden=true)
    % Defines the possible methods for FRET calculcation. The descriptions are
    % used when asking the user to specify what calculation to be used and what
    % assumptions can be made.
    % Hidden just to keep down the clutter when displaying objects.
    fretGeometryNames = {'tandem3','independent3','acceptor/total'};
    fretGeometryDesc  = {'Tandem 3-color (D->A1->A2)','Independent 3-color (D->A1, D->A2)','Acceptor/Total'};
end

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
     
    
     function changed = verifyFretGeometry(this)
     % Ensure that the method for calculating FRET is specified for 3/4-color
     % FRET. If not or if invalid, ask the user to specify.
        changed = false;
     
        if ~isChannel(this,'acceptor2'), return; end
     
        if ~isfield(this.fileMetadata,'fretGeometry') || ...
           ~ismember(this.fileMetadata.fretGeometry,this.fretGeometryNames),

            fretGeometry = 'acceptor/total';  %default

            % Handle legacy tag that specifies tandem3 calculation.
            if isfield(this.fileMetadata,'isTandem3') && this.fileMetadata.isTandem3,
                fretGeometry = 'tandem3';

            % Otherwise, we have to ask the user what to do.
            else
                [sel,ok] = listdlg('PromptString','What is the 3-color FRET geometry in this experiment?', ...
                               'SelectionMode','single','ListSize',[300 120], ...
                               'ListString',TracesFret4.fretGeometryDesc );
                if ok,
                    fretGeometry = TracesFret4.fretGeometryNames{sel};
                end
            end

            this.fileMetadata.fretGeometry = fretGeometry;
            changed = true;
        end
     end
     
     
    %% ================ DATA MANIPULATION ================ %%
    
    function this = recalculateFret( this, thresholds, indexes )
    % Recalculate FRET efficiencies (all fields) using fluorescence. 
    %    data.recalculateFret( THRESH );
    % 
    % Replaces the current "fret" properties with the newly calculated traces.
    % Thresholds of total intensity below which FRET is set to zero can be given
    % as the extra parameter (THRESH, Nx1 array, one per trace). Otherwise, they
    % calculated from the data as a few standard deviations above background.
    % FIXME: specify that THRESH is NxM, where M is the number of FRET channels.
    %
    % NOTE: There are three ways to calculate 3-color FRET:
    %   1) Tandem acceptors (D->A1->A2): assuming there is no donor->acceptor2
    %        FRET. 
    %   2) Independent acceptors (D->A1, D->A2): assuming there is no
    %        acceptor1->acceptor2 FRET.
    %   2) Unknown: if the above assumptions can't be made, it isn't possible
    %        to calculate FRET. Instead, calculate the fraction of acceptor
    %        intensity vs total (A1/total, A2/total). This is not FRET!
    %        It just tells you which acceptor is brighter.
    % 
    % This function will look in fileMetadata to find out which to use, and if
    % not known, acceptor/total is assumed. If the metadata is not available in
    % a .traces file, the function that loads it should ask the user.
    % 2-color FRET will be calculated as usual (see TracesFret).
    % 
    
        %--------------  Determine how to calculate FRET  --------------%
        
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
        
        % If the method for calculating FRET is not given or is not valid, ask
        % the user for the correct method.
        this.verifyFretGeometry();
        
        
        %------------------  Actually calculate FRET  ------------------%
        if nargin<2,
            thresholds = zeros(this.nTraces,1);  %fret1 only
        else
            if size(thresholds,2)>1,
                warning('Only fret1 thresholds will be used');
            end
        end
        
        % Get list of traces to correct; correct all if not specified.
        if nargin<3,
            indexes = 1:this.nTraces;
        end
        if islogical(indexes),  indexes = find(indexes);  end
        if isempty(indexes), return; end
        assert( isvector(indexes) );
        indexes = to_row(indexes);
        
        % Calculate FRET for each trace.
        assert( this.isChannel('fret') && this.isChannel('fret2') );
        total = this.total;
        
        switch this.fileMetadata.fretGeometry,
        case 'tandem3',
            % D1->A1->A2 FRET. Assumes no donor->acceptor2 FRET.
            acc = this.acceptor + this.acceptor2; %total acceptor intensity
            this.fret  = acc./total;
            this.fret2 = this.acceptor2./acc;
        case 'independent3',
            % D1->A1, D1->A2 FRET. Assumes no acceptor1->acceptor2 FRET.
            this.fret  = this.acceptor ./(this.acceptor +this.donor);
            this.fret2 = this.acceptor2./(this.acceptor2+this.donor);
        case 'independent4',
            % D1->A1, D2->A2 FRET. Assumes the two pairs are fully independent,
            % which means no D1->A2, D2->A1, A1->A2, etc FRET.
            this.fret  = this.acceptor ./(this.acceptor +this.donor);
            this.fret2 = this.acceptor2./(this.acceptor2+this.donor2);
        case 'acceptor/total',
            % Fraction of total intensity in each channel (Not FRET).
            this.fret  = this.acceptor./total;
            this.fret2 = this.acceptor2./total;
        end
        
        
        %--------------  Apply blinking/bleaching corrections  --------------%
        
        % Set FRET to zero when donor is dark (below total intensity threshold).
        lt = max(1, calcLifetime(total) );  %point where donor bleaches in each trace.
        constants = cascadeConstants;
        
        for i=indexes,
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
            if strcmp(this.fileMetadata.fretGeometry,'tandem3'),
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






