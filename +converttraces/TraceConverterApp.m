classdef TraceConverterApp < matlab.apps.AppBase
% Author: Zeliha Kilic
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % App Properties
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties (Access = public)
        UIFigure matlab.ui.Figure
        DeepLassiWarningShown logical = false;

        % Tab Group
        TabGroup matlab.ui.container.TabGroup
        
        % Tab 1: OpenFRET
        OpenFRETTab matlab.ui.container.Tab
        OpenFRETSelectButton matlab.ui.control.Button
        OpenFRETFileList matlab.ui.control.ListBox
        OpenFRETTitleField matlab.ui.control.EditField
        OpenFRETDescriptionField matlab.ui.control.EditField
        OpenFRETDateField matlab.ui.control.EditField
        OpenFRETExperimentIDField matlab.ui.control.EditField
        OpenFRETAuthorsField matlab.ui.control.EditField
        OpenFRETInstitutionField matlab.ui.control.EditField
        OpenFRETMicroscopeField matlab.ui.control.EditField
        OpenFRETDetectorField matlab.ui.control.EditField
        OpenFRETObjectiveField matlab.ui.control.EditField
        OpenFRETWavelength1Field matlab.ui.control.NumericEditField
        OpenFRETWavelength2Field matlab.ui.control.NumericEditField
        OpenFRETUseChannel1 matlab.ui.control.CheckBox
        OpenFRETUseChannel2 matlab.ui.control.CheckBox
        OpenFRETCrosstalkField matlab.ui.control.NumericEditField
        OpenFRETConvertButton matlab.ui.control.Button
        OpenFRETLogArea matlab.ui.control.TextArea
        
        % Tab 2: BNPhdf5
        BNPTab matlab.ui.container.Tab
        BNPSelectButton matlab.ui.control.Button
        BNPFileList matlab.ui.control.ListBox
        BNPConvertButton matlab.ui.control.Button
        BNPLogArea matlab.ui.control.TextArea
        
        % Tab 3: tmaven hdf5
        TmavenTab matlab.ui.container.Tab
        TmavenSelectButton matlab.ui.control.Button
        TmavenFileList matlab.ui.control.ListBox
        TmavenConvertButton matlab.ui.control.Button
        TmavenLogArea matlab.ui.control.TextArea
        
        % Tab 4: .mat
        MatTab matlab.ui.container.Tab
        MatSelectButton matlab.ui.control.Button
        MatFileList matlab.ui.control.ListBox
        MatConvertButton matlab.ui.control.Button
        MatLogArea matlab.ui.control.TextArea

        % Tab 5: 3 color Alex .mat, .npz
        DLassiTab matlab.ui.container.Tab
        DLassiSelectButton matlab.ui.control.Button
        DLassiFileList matlab.ui.control.ListBox
        DLassiFilePath matlab.ui.control.Button
        DLassiFilePathList matlab.ui.control.ListBox
        PyEnvPath matlab.ui.control.Button
        PyEnvPathFileList matlab.ui.control.ListBox
        DLassiConvertButton matlab.ui.control.Button
        DLassiLogArea matlab.ui.control.TextArea
        DLassiPathLogArea matlab.ui.control.TextArea
        PyEnvPathLogArea matlab.ui.control.TextArea

    end

    properties (Access = private)
        % File lists for each tab
        OpenFRETFiles cell
        BNPFiles cell
        TmavenFiles cell
        MatFiles cell
        DLassiFiles cell
        DLassiPathFiles cell
        PyEnvPathFiles cell
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % App Methods
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function onTabChanged(app, ~, event)
            if event.NewValue == app.DLassiTab && ~app.DeepLassiWarningShown
                uialert(app.UIFigure, ...
                    ['DeepLASI conversion requires DeepLASI already installed.' newline ...
                    'Please select a valid DeepLASI data/output folder and Python executable.'], ...
                    'DeepLASI Prerequisite', ...
                    'Icon', 'warning');
                app.DeepLassiWarningShown = true;
            end
        end
        %----------------------------------------------------------
        function log(app, logArea, msg)
            % Log message to specified text area
            timestamp = datestr(now, 'HH:MM:SS');
            logArea.Value{end+1} = sprintf('[%s] %s', timestamp, msg);
            drawnow;
        end

        %----------------------------------------------------------
        function SelectFiles(app, tabName)
            % Generic file selection for all tabs
            [f, p] = uigetfile( ...
                {'*.traces', 'SPARTAN .traces'; ...
                '*.rawtraces', 'SPARTAN Raw .rawtraces'; ...
                '*.*', 'All files'}, ...
                 'Select trace files', ...
                'MultiSelect', 'on');
            if isequal(f, 0)
                return;
            end
            
            if iscell(f)
                files = cellfun(@(x) fullfile(p, x), f, 'UniformOutput', false);
            else
                files = {fullfile(p, f)};
            end
            
            % Store files and update UI based on tab
            switch tabName
                case 'OpenFRET'
                    app.OpenFRETFiles = files;
                    app.OpenFRETFileList.Items = files;
                    app.log(app.OpenFRETLogArea, sprintf('Selected %d file(s)', numel(files)));
                case 'BNP'
                    app.BNPFiles = files;
                    app.BNPFileList.Items = files;
                    app.log(app.BNPLogArea, sprintf('Selected %d file(s)', numel(files)));
                case 'Tmaven'
                    app.TmavenFiles = files;
                    app.TmavenFileList.Items = files;
                    app.log(app.TmavenLogArea, sprintf('Selected %d file(s)', numel(files)));
                case 'Mat'
                    app.MatFiles = files;
                    app.MatFileList.Items = files;
                    app.log(app.MatLogArea, sprintf('Selected %d file(s)', numel(files)));
                case 'DLassi'
 
                    app.DLassiFiles = files;
                    app.DLassiFileList.Items = files;
                    app.log(app.DLassiLogArea, sprintf('Selected %d file(s)', numel(files)));                    
            
            
            end
        end
        %----------------------------------------------------------

        function SelectPyEnvPath(app, tabName)
            % Generic file selection for all tabs
            [f,p] = uigetfile('Select python env path');
            
            if iscell(f)
                files = cellfun(@(x) fullfile(p, x), f, 'UniformOutput', false);
            else
                files = {fullfile(p, f)};
            end
            
            % Store files and update UI based on tab
            switch tabName

                case 'DLassi'

                    app.PyEnvPathFiles =files;
                    app.log(app.DLassiLogArea, sprintf('Selected %d path(s)', numel(app.PyEnvPathFiles)));

            end
        end

        %----------------------------------------------------------

        function SelectDLassiFilePath(app, tabName)


            % Generic file selection for all tabs
            [p] = uigetdir('Select Deep Lasi data path');

            if iscell(p)
                files = cellfun(@(x) x, p, 'UniformOutput', false);
            else
                files = {p};
            end
            
            % Store files and update UI based on tab
            switch tabName

                case 'DLassi'
                    app.DLassiPathFiles =files;
                    app.log(app.DLassiLogArea, sprintf('Selected %d path(s)', numel(files)));

            end
        end

        %----------------------------------------------------------
        function ConvertToOpenFRET(app, ~)
            if isempty(app.OpenFRETFiles)
                uialert(app.UIFigure, 'Please select at least one .traces file first.', 'No Files Selected');
                return;
            end
            
            app.log(app.OpenFRETLogArea, 'Starting OpenFRET conversion...');
            
            try
                % Get metadata from UI
                title = app.OpenFRETTitleField.Value;
                description = app.OpenFRETDescriptionField.Value;
                date = app.OpenFRETDateField.Value;
                experiment_id = app.OpenFRETExperimentIDField.Value;
                authors_str = app.OpenFRETAuthorsField.Value;
                institution = app.OpenFRETInstitutionField.Value;
                microscope = app.OpenFRETMicroscopeField.Value;
                detector = app.OpenFRETDetectorField.Value;
                objective = app.OpenFRETObjectiveField.Value;
                
                % Parse authors (comma-separated)
                if ~isempty(authors_str)
                    authors = strsplit(authors_str, ',');
                    authors = cellfun(@strtrim, authors, 'UniformOutput', false);
                else
                    authors = {'User'};
                end
                
                excitation_wavelength.channel1 = app.OpenFRETWavelength1Field.Value;
                excitation_wavelength.channel2 = app.OpenFRETWavelength2Field.Value;
                
                options.use_channel1 = app.OpenFRETUseChannel1.Value;
                options.use_channel2 = app.OpenFRETUseChannel2.Value;
                options.donor_crosstalk = app.OpenFRETCrosstalkField.Value;
                
                % Process each file
                for i = 1:numel(app.OpenFRETFiles)
                    app.log(app.OpenFRETLogArea, sprintf('Processing file %d of %d: %s', i, numel(app.OpenFRETFiles), app.OpenFRETFiles{i}));
                    
                    inputFile = app.OpenFRETFiles{i};
                    [inputPath, inputName, ~] = fileparts(inputFile);
                    
                    % Check if input is on network drive (Z: or other network paths)
                    % Save to local directory first, then copy back to Z drive
                    % in linux machine unless you paste the networked
                    % copied to desktop, you cannot upload it to the
                    % website,
                    % you can upload the Z dirve file form windows though.
                    % To make sure all is well, I just have two ways of
                    % saving the files.
                    if startsWith(inputPath, 'Z:') || startsWith(inputPath, '\\') || ...
                       (isunix && ~startsWith(inputPath, '/home') && ~startsWith(inputPath, '/tmp'))
                        % Use temp directory for local save
                        localDir = fullfile(tempdir, 'OpenFRET_conversions');
                        if ~exist(localDir, 'dir')
                            mkdir(localDir);
                        end
                        
                        app.log(app.OpenFRETLogArea, sprintf('Network drive detected. Using local temp: %s', localDir));
                        
                        % Create local copy of input file
                        localInputFile = fullfile(localDir, [inputName, '.traces']);
                        app.log(app.OpenFRETLogArea, sprintf('Step 1: Copying input to local: %s', localInputFile));
                        copyfile(inputFile, localInputFile,'f');
                        
                        % Convert using local file (output will be saved in local directory)
                        app.log(app.OpenFRETLogArea, sprintf('Step 2: Converting in local directory...'));
                        converttraces.convertToOpenFRET({localInputFile}, ...
                            'title', title, ...
                            'description', description, ...
                            'date', date, ...
                            'experiment_id', experiment_id, ...
                            'authors', authors, ...
                            'institution', institution, ...
                            'microscope', microscope, ...
                            'detector', detector, ...
                            'objective', objective, ...
                            'excitation_wavelength', excitation_wavelength, ...
                            'options', options);
                        
                        % Output files are saved in local directory - leave them there
                        localZipFile = fullfile(localDir, [inputName, '.json.zip']);

                        % Here I copy them back to the original location. I found out that whne I create the files in the network drive, then they dont become compatible with the server
                        % Original input directory (this is teh file submitted location, not the local temp)
                        
                        inputOpenFRETDir = fullfile(inputPath, 'OpenFRET_conversions');
                        if ~exist(inputOpenFRETDir, 'dir')
                            mkdir(inputOpenFRETDir);
                        end
                            
                        inputDir = inputOpenFRETDir; % Save in OpenFRET_conversions subfolder on network drive
                        % Destination paths
                        destZipFile  = fullfile(inputDir, [inputName, '.json.zip']);


                        if isfile(localZipFile)
                            copyfile(localZipFile, destZipFile, 'f');
                        end
                        

                        if exist(localZipFile, 'file')
                            app.log(app.OpenFRETLogArea, sprintf('Zip file saved to local directory: %s', localZipFile));
                        end
                        
                        app.log(app.OpenFRETLogArea, sprintf('Files saved in local directory: %s', localDir));
                        
                        % Clean up only the input file copy, keep output files
                        if exist(localInputFile, 'file'), delete(localInputFile); end
                        % Now I can get rid of the local copy of that file
                        delete(localZipFile)
                        
                    else
                        % File is already on local drive, convert directly
                        converttraces.convertToOpenFRET({inputFile}, ...
                            'title', title, ...
                            'description', description, ...
                            'date', date, ...
                            'experiment_id', experiment_id, ...
                            'authors', authors, ...
                            'institution', institution, ...
                            'microscope', microscope, ...
                            'detector', detector, ...
                            'objective', objective, ...
                            'excitation_wavelength', excitation_wavelength, ...
                            'options', options);
                    end
                end
                
                app.log(app.OpenFRETLogArea, sprintf('Successfully converted %d file(s) to OpenFRET format.', numel(app.OpenFRETFiles)));
                uialert(app.UIFigure, 'Conversion completed successfully!', 'Success', 'Icon', 'success');
                
            catch ME
                app.log(app.OpenFRETLogArea, sprintf('ERROR: %s', ME.message));
                uialert(app.UIFigure, sprintf('Conversion failed: %s', ME.message), 'Error', 'Icon', 'error');
            end
        end

        %----------------------------------------------------------
        function ConvertToBNP(app, ~)
            if isempty(app.BNPFiles)
                uialert(app.UIFigure, 'Please select at least one .traces file first.', 'No Files Selected');
                return;
            end
            
            app.log(app.BNPLogArea, 'Starting BNP hdf5 conversion...');
            
            try
                for i = 1:numel(app.BNPFiles)
                    app.log(app.BNPLogArea, sprintf('Processing file %d of %d: %s', i, numel(app.BNPFiles), app.BNPFiles{i}));
                    converttraces.TracesToBNPh5(app.BNPFiles{i});
                end
                
                app.log(app.BNPLogArea, sprintf('Successfully converted %d file(s) to BNP hdf5 format.', numel(app.BNPFiles)));
                uialert(app.UIFigure, 'Conversion completed successfully!', 'Success', 'Icon', 'success');
                
            catch ME
                app.log(app.BNPLogArea, sprintf('ERROR: %s', ME.message));
                uialert(app.UIFigure, sprintf('Conversion failed: %s', ME.message), 'Error', 'Icon', 'error');
            end
        end

        %----------------------------------------------------------
        function ConvertToTmaven(app, ~)
            if isempty(app.TmavenFiles)
                uialert(app.UIFigure, 'Please select at least one .traces file first.', 'No Files Selected');
                return;
            end
            
            app.log(app.TmavenLogArea, 'Starting tmaven hdf5 conversion...');
            
            try
                for i = 1:numel(app.TmavenFiles)
                    app.log(app.TmavenLogArea, sprintf('Processing file %d of %d: %s', i, numel(app.TmavenFiles), app.TmavenFiles{i}));
                    converttraces.TracesTohdf5(app.TmavenFiles{i});
                end
                
                app.log(app.TmavenLogArea, sprintf('Successfully converted %d file(s) to tmaven hdf5 format.', numel(app.TmavenFiles)));
                uialert(app.UIFigure, 'Conversion completed successfully!', 'Success', 'Icon', 'success');
                
            catch ME
                app.log(app.TmavenLogArea, sprintf('ERROR: %s', ME.message));
                uialert(app.UIFigure, sprintf('Conversion failed: %s', ME.message), 'Error', 'Icon', 'error');
            end
        end

        %----------------------------------------------------------
        function ConvertToMat(app, ~)
            if isempty(app.MatFiles)
                uialert(app.UIFigure, 'Please select at least one .traces file first.', 'No Files Selected');
                return;
            end
            
            app.log(app.MatLogArea, 'Starting .mat conversion...');
            
            try
                for i = 1:numel(app.MatFiles)
                    app.log(app.MatLogArea, sprintf('Processing file %d of %d: %s', i, numel(app.MatFiles), app.MatFiles{i}));
                     converttraces.tracesToMat(app.MatFiles{i});
                end
                
                app.log(app.MatLogArea, sprintf('Successfully converted %d file(s) to .mat format.', numel(app.MatFiles)));
                uialert(app.UIFigure, 'Conversion completed successfully!', 'Success', 'Icon', 'success');
                
            catch ME
                app.log(app.MatLogArea, sprintf('ERROR: %s', ME.message));
                uialert(app.UIFigure, sprintf('Conversion failed: %s', ME.message), 'Error', 'Icon', 'error');
            end
        end
        %----------------------------------------------------------
        function ConvertToDLassiAlex3Color(app, ~)


            if isempty(app.DLassiFiles)
                uialert(app.UIFigure, 'Please select at least one .traces file first.', 'No Files Selected');
                return;
            end
            app.log(app.DLassiLogArea, 'Starting Deep Lasi npz,mat conversion...');

            if isempty(app.PyEnvPath)
                uialert(app.UIFigure, 'Please select the Python eneviroment file first.', 'No Python Path Selected');
                return;
            end

            if isempty(app.DLassiFilePath)
                uialert(app.UIFigure, 'Please select the Deep Lasi data path first.', 'No Deep LAsi /data Path Selected');
                return;
            end
                        
            try
                for i = 1:numel(app.DLassiFiles)
                    app.log(app.DLassiLogArea, sprintf('Processing file %d of %d: %s', i, numel(app.DLassiFiles), app.DLassiFiles{i}));
                    converttraces.tracesToNpz(app.DLassiFiles{i},app.DLassiPathFiles{1},app.PyEnvPathFiles{1});
                end
                
                app.log(app.DLassiLogArea, sprintf('Successfully converted %d file(s) to Deep Lasi .npz,.mat format.', numel(app.DLassiFiles)));
                uialert(app.UIFigure, 'Conversion completed successfully!', 'Success', 'Icon', 'success');
                
            catch ME
                app.log(app.DLassiLogArea, sprintf('ERROR: %s', ME.message));
                uialert(app.UIFigure, sprintf('Conversion failed: %s', ME.message), 'Error', 'Icon', 'error');
            end
        end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % UI Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = private)

        function createComponents(app)
            % Create main figure - wider for single column layout
            app.UIFigure = uifigure('Position', [100 100 710 500], ...
                'Name', 'Trace Converter', ...
                'Resize', 'on');
            
            % Create tab group
            app.TabGroup = uitabgroup(app.UIFigure, 'Position', [10 10 690 480]);
            
            % Initialize file lists
            app.OpenFRETFiles = {};
            app.BNPFiles = {};
            app.TmavenFiles = {};
            app.MatFiles = {};
            app.DLassiFiles = {};
            app.DLassiPathFiles = {};       
            app.PyEnvPathFiles = {};     

            % Create tabs
            createOpenFRETTab(app);
            createBNPTab(app);
            createTmavenTab(app);
            createMatTab(app);
            createNpzTab(app);
        end

        function createOpenFRETTab(app)
            % Tab 1: OpenFRET - Single column layout
            app.OpenFRETTab = uitab(app.TabGroup, 'Title', '.traces to OpenFRET');
            
            % LAYOUT CONFIGURATION
            % Top section
            topY = 420;                    % Y position for select button
            buttonWidth = 200;             % Width of buttons
            buttonHeight = 30;             % Height of buttons
            tabWidth = 530;                % Total tab width
            buttonX = (tabWidth - buttonWidth) / 2;  % Center button horizontally
            
            % File list section
            fileListLabelY = 390;          % Y position for "Selected Files" label
            fileListY = 350;               % Y position for file list
            fileListHeight = 30;           % Height of file list
            fileListWidth = 490;           % Width of file list
            leftMargin = 20;               % Left margin for labels/fields
            
            % Metadata fields section
            yPos = 320;                    % Starting Y position for metadata fields
            spacing = 25;                  % Vertical spacing between fields
            labelWidth = 120;              % Width of labels
            fieldX = 150;                  % X position for input fields
            fieldWidth = 410;              % Width of input fields
            
            % Bottom section
            logLabelY = 50;                % Y position for "Log:" label
            logY = 20;                     % Y position for log area
            logHeight = 30;                % Height of log area
            logWidth = 500;                % Width of log area
            % 
            
            % File selection button (centered)
            app.OpenFRETSelectButton = uibutton(app.OpenFRETTab, ...
                'Text', 'Select .traces Files', ...
                'Position', [buttonX topY buttonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.SelectFiles('OpenFRET'));
            
            % Small window to see selected traces
            uilabel(app.OpenFRETTab, 'Text', 'Selected Files:', ...
                'Position', [leftMargin fileListLabelY 100 22]);
            app.OpenFRETFileList = uilistbox(app.OpenFRETTab, ...
                'Position', [leftMargin fileListY fileListWidth fileListHeight], ...
                'Items', {});
            
            % Metadata fields
            uilabel(app.OpenFRETTab, 'Text', 'Title:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETTitleField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', '(Experiment name goes here)');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Description:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETDescriptionField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', '(Short description of instrument goes here)');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Date:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETDateField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos 200 22], 'Value', 'YYYY-MM-DD');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Experiment ID:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETExperimentIDField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos 200 22], 'Value', 'YYYYMMDD_ABC_1');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Authors (comma-separated):', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETAuthorsField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', '(Author of experiment)');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Institution:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETInstitutionField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', 'St. Jude Childrens Research Hospital');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Microscope:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETMicroscopeField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', 'TIRF 2/3');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Detector:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETDetectorField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', 'Kinetix scMOS');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Objective:', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETObjectiveField = uieditfield(app.OpenFRETTab, 'text', ...
                'Position', [fieldX yPos fieldWidth 22], 'Value', '60X 1.5 NA');
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Excitation Wavelength 1 (nm):', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETWavelength1Field = uieditfield(app.OpenFRETTab, 'numeric', ...
                'Position', [fieldX yPos 100 22], 'Value', 532);
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Excitation Wavelength 2 (nm):', 'Position', [leftMargin yPos labelWidth 22]);
            app.OpenFRETWavelength2Field = uieditfield(app.OpenFRETTab, 'numeric', ...
                'Position', [fieldX yPos 100 22], 'Value', 640);
            yPos = yPos - spacing;
            
            % Options
            app.OpenFRETUseChannel1 = uicheckbox(app.OpenFRETTab, ...
                'Text', 'Use Channel 1 (Donor)', ...
                'Position', [leftMargin+300 yPos+50 200 22], 'Value', true);
            yPos = yPos - spacing;
            
            app.OpenFRETUseChannel2 = uicheckbox(app.OpenFRETTab, ...
                'Text', 'Use Channel 2 (Acceptor)', ...
                'Position', [leftMargin+300 yPos+50 200 22], 'Value', true);
            yPos = yPos - spacing;
            
            uilabel(app.OpenFRETTab, 'Text', 'Donor Crosstalk:', 'Position', [leftMargin yPos+50 labelWidth 22]);
            app.OpenFRETCrosstalkField = uieditfield(app.OpenFRETTab, 'numeric', ...
                'Position', [fieldX yPos+50 100 22], 'Value', 0.09);
            yPos = yPos - spacing;
            
            % Convert button (centered)
            convertButtonWidth = 130;
            convertButtonX =( (tabWidth - convertButtonWidth) / 2 )+180;
            app.OpenFRETConvertButton = uibutton(app.OpenFRETTab, ...
                'Text', 'Convert to OpenFRET', ...
                'Position', [convertButtonX-50 yPos+70 convertButtonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.ConvertToOpenFRET,'BackgroundColor',[209 25 71]/255);
            
            % Small log area at bottom
            uilabel(app.OpenFRETTab, 'Text', 'Log:', 'Position', [leftMargin logLabelY-20 50 22]);
            app.OpenFRETLogArea = uitextarea(app.OpenFRETTab, ...
                'Position', [leftMargin logY-15 logWidth logHeight-5], ...
                'Editable', 'off');
        end

        function createBNPTab(app)
            % Tab 2: BNP hdf5 - Single column layout
            app.BNPTab = uitab(app.TabGroup, 'Title', '.traces to BNPhdf5');
            
            tabWidth = 530;
            convertButtonWidth = 150;
            convertButtonX =( (tabWidth - convertButtonWidth) / 2 )+180;
            buttonHeight = 30;             % Height of buttons

            % File selection button
            app.BNPSelectButton = uibutton(app.BNPTab, ...
                'Text', 'Select .traces Files', ...
                'Position', [165 420 200 30], ...
                'ButtonPushedFcn', @(~,~) app.SelectFiles('BNP'));
            
            % Small window to see selected traces (moved down)
            uilabel(app.BNPTab, 'Text', 'Selected Files:', ...
                'Position', [20 390 100 22]);
            app.BNPFileList = uilistbox(app.BNPTab, ...
                'Position', [20 340 490 50], ...
                'Items', {});
            
            % Convert button (centered)
            app.BNPConvertButton = uibutton(app.BNPTab, ...
                'Text', 'Convert to BNPhdf5', ...
                'Position', [convertButtonX 100 convertButtonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.ConvertToBNP,'BackgroundColor',[209 25 71]/255);
            
            % Small log area at bottom
            uilabel(app.BNPTab, 'Text', 'Log:', 'Position', [20 80-30 50 22]);
            app.BNPLogArea = uitextarea(app.BNPTab, ...
                'Position', [20 20 490 30], ...
                'Editable', 'off');
        end

        function createTmavenTab(app)
            % Tab 3: tmaven hdf5 - Single column layout
            app.TmavenTab = uitab(app.TabGroup, 'Title', '.traces to tmaven hdf5');
             tabWidth = 530;
            convertButtonWidth = 150;
            convertButtonX =( (tabWidth - convertButtonWidth) / 2 )+180;
            buttonHeight = 30;    
            % File selection button (centered, moved down)
            app.TmavenSelectButton = uibutton(app.TmavenTab, ...
                'Text', 'Select .traces Files', ...
                'Position', [165 420 200 30], ...
                'ButtonPushedFcn', @(~,~) app.SelectFiles('Tmaven'));
            
            % Small window to see selected traces (moved down)
            uilabel(app.TmavenTab, 'Text', 'Selected Files:', ...
                'Position', [20 390 100 22]);
            app.TmavenFileList = uilistbox(app.TmavenTab, ...
                'Position', [20 340 490 50], ...
                'Items', {});
            
            % Convert button (centered)
            app.TmavenConvertButton = uibutton(app.TmavenTab, ...
                'Text', 'Convert to tmaven hdf5', ...
                'Position', [convertButtonX 100 convertButtonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.ConvertToTmaven,'BackgroundColor',[209 25 71]/255);
            
            % Small log area at bottom
            uilabel(app.TmavenTab, 'Text', 'Log:', 'Position', [20 80-30 50 22]);
            app.TmavenLogArea = uitextarea(app.TmavenTab, ...
                'Position', [20 20 490 30], ...
                'Editable', 'off');
        end

        function createMatTab(app)
            % Tab 4: .mat - Single column layout
            app.MatTab = uitab(app.TabGroup, 'Title', '.traces to .mat');
            tabWidth = 530;
            convertButtonWidth = 150;
            convertButtonX =( (tabWidth - convertButtonWidth) / 2 )+180;
            buttonHeight = 30;    
            % File selection button (centered, moved down)
            app.MatSelectButton = uibutton(app.MatTab, ...
                'Text', 'Select .traces Files', ...
                'Position', [165 420 200 30], ...
                'ButtonPushedFcn', @(~,~) app.SelectFiles('Mat'));
            
            % Small window to see selected traces (moved down)
            uilabel(app.MatTab, 'Text', 'Selected Files:', ...
                'Position', [20 390 100 22]);
            app.MatFileList = uilistbox(app.MatTab, ...
                'Position', [20 340 490 50], ...
                'Items', {});
            
            % Convert button (centered)
            app.MatConvertButton = uibutton(app.MatTab, ...
                'Text', 'Convert to .mat', ...
                'Position', [convertButtonX 100 convertButtonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.ConvertToMat,'BackgroundColor',[209 25 71]/255);
            
            % Small log area at bottom
            uilabel(app.MatTab, 'Text', 'Log:', 'Position', [20 80-30 50 22]);
            app.MatLogArea = uitextarea(app.MatTab, ...
                'Position', [20 20 490 30], ...
                'Editable', 'off');
        end




        function createNpzTab(app)

        % after creating TabGroup
        app.TabGroup.SelectionChangedFcn = @(src,event) app.onTabChanged(src,event);

            % Tab 2: Deep LAssi Npz Tab - Single column layout
            app.DLassiTab = uitab(app.TabGroup, 'Title', '.traces to DeepLasiNpz');
            
            tabWidth = 530;
            convertButtonWidth = 150;
            convertButtonX =( (tabWidth - convertButtonWidth) / 2 )+180;
            buttonHeight = 20;             % Height of buttons
            shifty = -50;
            shifty1 = -50;


            % File selection button


            app.DLassiSelectButton = uibutton(app.DLassiTab, ...
                'Text', 'Select .traces/.rawtraces Files', ...
                'Position', [165 420+shifty 250 20], ...
                'ButtonPushedFcn', @(~,~) app.SelectFiles('DLassi'));

            app.DLassiFilePath = uibutton(app.DLassiTab, ...
                'Text', 'Select DeepLASI Simulation Data Path ', ...
                'Position', [165 450+shifty 250 20], ...
                'ButtonPushedFcn', @(~,~) app.SelectDLassiFilePath('DLassi')); 

            app.PyEnvPath = uibutton(app.DLassiTab, ...
                'Text', 'Select Python Env Path ', ...
                'Position', [165 480+shifty 250 20], ...
                'ButtonPushedFcn', @(~,~) app.SelectPyEnvPath('DLassi')); 


            % Small window to see selected traces (moved down)
            uilabel(app.DLassiTab, 'Text', 'Selected Files:', ...
                'Position', [20 370 100 22]);
            app.DLassiFileList = uilistbox(app.DLassiTab, ...
                'Position', [20 320+shifty1 490 30], ...
                'Items', {});
   

            % Convert button (centered)
            app.DLassiConvertButton = uibutton(app.DLassiTab, ...
                'Text', 'Convert to Deep Lasi npz', ...
                'Position', [convertButtonX 180 convertButtonWidth buttonHeight], ...
                'ButtonPushedFcn', @(~,~) app.ConvertToDLassiAlex3Color,'BackgroundColor',[209 25 71]/255);
            
            % Small log area at bottom
            uilabel(app.DLassiTab, 'Text', 'Log:', 'Position', [20 200-30 50 22]);
            app.DLassiLogArea = uitextarea(app.DLassiTab, ...
                'Position', [20 110 500 60], ...
                'Editable', 'off');
  
        end

    end



   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % App Constructor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Access = public)

        function app = TraceConverterApp
            createComponents(app);
        end
    end
end

