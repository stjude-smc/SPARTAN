function USERParams = MetaDataOpenFRET 

USERParams.TITLE = "(Experiment name goes here)";
USERParams.DESCRIPTION = "(Short description of instrument goes here)";
USERParams.EXPERIMENT_TYPE = "(Experiment date goes here)";
USERParams.AUTHORS = "(Author of experiment)";
USERParams.INSTITUTION = "St. Jude Chidren's Research HOSPITAL ";
USERParams.DATE = "YYYY-MM-DD";
USERParams.EXPERIMENT_ID = "YYYYMMDD_ABC_1"; % Recommended format: year_initials_experiment
USERParams.BUFFER_CONDITIONS = "4X PBS";
USERParams.TEMPERATURE = "28 C";
USERParams.MICROSCOPE = "TIRF 2/3";
USERParams.DETECTOR = "KINETIX";
USERParams.OBJECTIVE = "60X 1.5 NA";
USERParams.EXCITATION_WAVELENGTH.CHANNEL1 = "532"; % In nanometers
USERParams.EXCITATION_WAVELENGTH.CHANNEL2 = "640"; % In nanometers

USERParams.OPTIONS.USE_CHANNEL1 = true; % true = include channel 1 (green/donor) in OpenFRET; false = omit
USERParams.OPTIONS.USE_CHANNEL2  = true; % true = include channel 2 (red/acceptor) in OpenFRET; false = omit

USERParams.OPTIONS.DONOR_CROSSTALK = 0.09; % (For 2-channel traces files): crosstalk of donor into acceptor channel, as fraction of donor-channel signal

end