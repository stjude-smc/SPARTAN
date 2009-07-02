function accTime_script(filenames, titles)




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
makeplotsOptions.pophistSumlen = 600; %2 min. at 200ms
% makeplotsOptions.cplot_scale_factor = 15; %2 min. at 200ms

selectionCriteria.min_lifetime  = 300; %donor lifetime
selectionCriteria.min_snr       = 8; 
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
    makeplots( autoFilenames(i), titles(i), 'targetAxes', {h} );
    xlim([0 makeplotsOptions.pophistSumlen]);
    set(gca, 'xtick', [0:150:600] );
    set(gca, 'xticklabel',{'0','30','60','90','120'} );
    xlabel('Time (sec)');
    if i>1, ylabel(''); end
    %ylim([-0.1 0.8]);
end


%%
subplot(1,nFiles+1,nFiles+1);
accTime( autoFilenames, titles );
ylim([0 0.8]);



