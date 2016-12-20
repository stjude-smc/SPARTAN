function defaultFigLayout(hFig, cbNew, cbOpen, fhExport, editCallbacks)
% Add menus and toolbars for standard SPARTAN figures.
%
%   defaultFigLayout(FIG) modifies the appearance of the figure FIG to match
%   the style of typical SPARTAN functions. Some rarely used menus and toolbar
%   buttons are removed.
%
%   defaultFigLayout(FIG, MFILE) changes the callback functions for the New
%   and Open toolbar and menu options to call the function specified in
%   the string MFILE. The target function should take a figure handle as its
%   first argument.
%
%   defaultFigLayout(..., fhExport) adds an "Export as text" menu option that
%   calls the given function handle.
%
%   defaultFigLayout(..., fhExport, editCallbacks) adds menu options in the
%   Edit menu. Cell columns are menu text and function handle for callback.
%
%   This function is internal. Its behavior may change in future versions.
%
%   See also makeplots, occtime, avgFretTime.

%   Copyright 2016 Cornell University All Rights Reserved.

% Could be used without the fig argument to creat the new figure.
% mfile could be detected from dbstack automatically. Could also be a function
% handle.

% NOTE: could have some settings here for global customization, like removing
% all figure editing tools for simplicity, setable in cascadeConstants.


narginchk(1,5);


%% Hide some rarely-used menus/buttons to reduce clutter.
toHide = {'figMenuDesktop','figMenuWindow','DataManager.Linking','figMenuGenerateCode'...
    'figMenuFilePreferences','figMenuFileSaveWorkspaceAs','figMenuFileImportData', ...
    'figMenuEditFindFiles','figMenuEditClearWorkspace','figMenuEditClearCmdHistory', ...
    'figMenuEditClearCmdWindow','Plottools.PlottoolsOn','Plottools.PlottoolsOff'};
cellfun( @(x)set(findall(hFig,'tag',x),'Visible','off'), toHide );



%% Customize New and Open tools to call the target function instead.
if nargin>=2 && ~isempty(cbNew),
    hMenu = findall(hFig,'tag','figMenuUpdateFileNew');
    delete(allchild(hMenu));  %Remove New Script, Figure, etc options.
    set(hMenu, 'Callback',cbNew, 'Visible','on');
    hMenu = findall(hFig,'tag','Standard.NewFigure');
    set(hMenu, 'ClickedCallback',cbNew, 'Visible','on');
end

if nargin>=3 && ~isempty(cbOpen),
    hMenu = findall(hFig,'tag','Standard.FileOpen');
    set(hMenu, 'ClickedCallback',cbOpen );
    set( findall(hFig,'tag','figMenuOpen'), 'Callback',cbOpen );
end

if nargin>=4 && ~isempty(fhExport),
    hMenu = findall(hFig,'tag','figMenuGenerateCode');
    set(hMenu, 'Label','Export .txt files', 'Callback',fhExport, 'Visible','on');
end

% Add custom Edit menu options
if nargin>=5 && ~isempty(editCallbacks),
    % Remove any old custom menu items.
    delete( findall(hFig,'tag','spartanMenuEdit') );
    
    hEditMenu = findall(hFig, 'tag','figMenuEdit');
    set(hEditMenu,'Visible','on');

    ncb = size(editCallbacks,1);
    for i=1:ncb,
        uimenu( hEditMenu, 'Label',editCallbacks{i,1}, 'Callback',editCallbacks{i,2}, ...
                'Tag','spartanMenuEdit' );
    end

    % Put the new ones at the top.
    c = allchild(hEditMenu);
    set(c(end),'Separator','on');
    set(hEditMenu,'Children',c([(1+ncb):end 1:ncb]));
end


end