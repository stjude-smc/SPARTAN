function accPlots(filenames, titles)




%%

if nargin<1,
    % Get filenames
    filenames = getFiles('*.traces');
end

autoFilenames = strrep(filenames,'.traces','_auto.txt');
nFiles = numel(filenames);

if nFiles<1,
    return;
end

figure;

%% Generate titles for figures
if nargin<2,
    titles = strrep(filenames,'_',' ');
    for i=1:nFiles,
        [p,titles{i}] = fileparts( titles{i} );
    end
end

% Setup settings
makeplotsOptions.constants = cascadeConstants;
makeplotsOptions.constants.cplot_scale_factor = 15;
makeplotsOptions.pophist_sumlen = 600; %2 min. at 200ms
makeplotsOptions.contour_bin_size  = 0.03;
makeplotsOptions.fretRange = [-0.1 0.8];

selectionCriteria.min_lifetime  = 600; %donor lifetime
selectionCriteria.min_snr       = 8; %SNR1
selectionCriteria.max_ncross    = 4; %Cy3 blinks
selectionCriteria.eq_overlap    = 0; %single molecules only
selectionCriteria.maxTotalSigma = 2; %total intensity w/i 2 sigma


for i=1:nFiles,
    
    % Filter traces file
    [d,a,f,ids,time] = loadTraces( filenames{i} );
    stats = traceStat(d,a,f);
    idxPicked = pickTraces( stats, selectionCriteria );
    d = d(idxPicked,:);
    a = a(idxPicked,:);
    f = f(idxPicked,:);
    ids = ids(idxPicked);
    
    saveTraces( autoFilenames{i}, 'txt', d,a,f,ids,time );
    
    % Display contour plots with extended axes
    h = subplot(1,nFiles+1,i);
    makeplotsOptions.targetAxes = {h};
    
    makeplots( autoFilenames(i), titles(i), makeplotsOptions );
    set(gca, 'xtick', [0:150:600] );
    set(gca, 'xticklabel',{'0','30','60','90','120'} );
    xlabel('Time (sec)');
    if i>1, ylabel(''); end
end


%%
subplot(1,nFiles+1,nFiles+1);
accTime( autoFilenames, titles );
ylim([0 0.8]);



