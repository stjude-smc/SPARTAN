function output = hmmopt(methodName, openDialog)
% Set default parameters values for various kinetic modeling methods.
% Yes we know the spelling of 'seperately' is not correct.
% We keep everything here to encourage consistency

% FIXME: also return fields and descriptions for settingdlg?
% Should this function also call settingdlg??
% Also consider handling simulation methods here as well...


narginchk(1,2);
nargoutchk(0,1);
assert(ischar(methodName), 'First input must be the method name');


%% Set default parameter values on first run
persistent allopt;

if isempty(allopt)
    % Segmental k-means (SKM)
    allopt.skm.maxItr     = 100;
    allopt.skm.quiet      = false;
    allopt.skm.convLL     = 1e-4;
    %allopt.skm.convGrad  = 1e-4;
    allopt.skm.fixRates   = false;
    allopt.skm.zeroEnd    = false;
    allopt.skm.seperately = true;

    % Baum-welch
    allopt.bw.maxItr      = 200;
    allopt.bw.verbose     = true;
    allopt.bw.convLL      = 1e-5;
    allopt.bw.convGrad    = 1e-5;
    allopt.bw.fixRates    = false;
    allopt.bw.seperately  = false;

    % ebFRET
    %allopt.eb.minStates   = 1;
    %allopt.eb.maxStates   = 5;
    allopt.eb.maxRestarts  = 2; %10;
    allopt.eb.threshold    = 1e-5;

    % MIL (Together)
    allopt.mil.maxItr      = 150;
    allopt.mil.verbose     = true;
    allopt.mil.convLL      = 1e-5;
    allopt.mil.convGrad    = 1e-5;
    allopt.mil.verbose     = true;
    allopt.mil.removeDarkState = true;
    allopt.mil.UseParallel = cascadeConstants('enable_parfor');
    
    % MIL (Separately)
    allopt.mil_sep.verbose = false;

    % MPL
    allopt.mpl.maxItr   = 200;
    allopt.mpl.verbose  = true;
    allopt.mpl.convLL   = 1e-6;
    allopt.mpl.convGrad = 1e-6;

    % HMJP
    allopt.hmjp.maxItr   = 1000;
    allopt.hmjp.alpha    = 2;     % Uniformization; determines further refinements of jump times within a frame period
    allopt.hmjp.beta     = 10;    % higher number = slow rates. 1/(beta*eta)=peak escape rate in prior
    allopt.hmjp.eta      = 2;     % gamma distribution shape parameter: 4=peaked prior, 2=exp prior.
    allopt.hmjp.HMC_eps  = 0.01;  % Hamiltonian Monte Carlo integration step size
    allopt.hmjp.HMC_L    = 50;    % Hamiltonian Monte Carlo number of Leap-frog integration steps.
end


%% Prompt user to change values, if requested
switch methodName
    case {'skm','Segmental k-Means'}
        name = 'skm';
        prompts = {'Max iterations','Analyze traces individually','LL Convergence'};
        fields  = {'maxItr','seperately','convLL'};

    case {'BWoptimize','Baum-Welch'}
        name = 'bw';
        prompts = {'Max iterations','LL Convergence','Grad. Convergence'};
        fields  = {'maxItr',        'convLL',        'convGrad'};

    case {'runEbFret','ebFRET'}
        name = 'eb';
        prompts = {'Max restarts','Convergence threshold'};   %'Min states','Max states',
        fields  = {'maxRestarts', 'threshold'};   %'minStates', 'maxStates',

    case {'milOptimize','MIL (Together)'}
        name = 'mil';
        prompts = {'Remove dark state (class 1)','Max iterations','LL Convergence','Grad. Convergence'};
        fields  = {'removeDarkState',  'maxItr',        'convLL',        'convGrad'};

    case {'milOptimizeSeparately','MIL (Separately)'}
        name = 'mil_sep';
        [prompts,fields] = deal({});

    case {'mplOptimize','MPL'}
        name = 'mpl';
        prompts = {'Max iterations','LL Convergence','Grad. Convergence'};
        fields  = {'maxItr',        'convLL',        'convGrad'};

    case {'runHMJP','HMJP'}
        name = 'hmjp';
        prompts = {'Alpha (uniformization)','Beta (rate prior distribution scale parameter)',...
                  'Eta (rate prior distribution shape parameter)','HMC integration step size',...
                  'HMC leap-frog steps','Max iterations'};
        fields  = {'alpha','beta','eta','HMC_eps','HMC_L','maxItr'};

    otherwise
        error('Invalid method name');
end

output = allopt.(name);


%% Prompt user to change values, if requested
if nargin >= 2 && numel(fields)>0 && openDialog 
    %output = mergestruct(output, optIn);
    output = settingdlg(output, fields, prompts);
    if ~isempty(output)
        allopt.(name) = output;
    end
end


end  %function





