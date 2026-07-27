path_home = 'X:\2024\work\MossyFibres'; % path to the main code folder
data_home = 'X:\MFB\MFB_AH_2023\Data'; % folder path to all datasets
savepath = 'X:\MFB'; % path to save results
addpath([path_home '/Init']);
addpath([path_home '/preprocessing']);
addpath([path_home '/preprocessing/Axon grouping']);
addpath([path_home '/preprocessing/Correlation data']);
addpath([path_home '/preprocessing/GC modeling']);
addpath([path_home '/Analysis']);
addpath([path_home '/Analysis/Utilities']);
addpath([path_home '/Analysis/Analysis_pipeline']);
addpath([path_home '/Analysis/Analysis_pipeline/Utilities']);

ColorCodes

files = { 
    '171212_16_19_37', 'Animal1'; ...
    '200130_13_21_13 FunctAcq', 'Animal2'; ...
    '200130_13_36_14 FunctAcq', 'Animal2'; ...
    '200130_13_49_09 FunctAcq', 'Animal2'; ...
    '200130_14_02_12 FunctAcq', 'Animal2'; ...
    '200130_14_15_24 FunctAcq', 'Animal2'; ...
    '200130_14_29_30 FunctAcq', 'Animal2'; ...
    '191018_13_39_41', 'Animal3'; ...
    '191018_13_56_55', 'Animal3'; ...
    '191018_14_30_00', 'Animal3'; ...
    '191018_14_11_33', 'Animal3'; ...
    '191209_13_44_12', 'Animal4'; ...
    '191209_14_04_14', 'Animal4'; ...
    '191209_14_18_13', 'Animal4'; ...
    '191209_14_32_39', 'Animal4'; ...
    '191209_14_46_58', 'Animal4'; ...
    '191209_15_01_22', 'Animal4' ...
};


A1_data = {
    '171212_16_19_37';
    };

A2_data = {
    '200130_13_21_13 FunctAcq';
    '200130_13_36_14 FunctAcq';
    '200130_13_49_09 FunctAcq';
    '200130_14_02_12 FunctAcq';
    '200130_14_15_24 FunctAcq';
    '200130_14_29_30 FunctAcq';
    };

A3_data = {
    '191018_13_39_41';
    '191018_13_56_55';
    '191018_14_30_00';
    '191018_14_11_33';
    };

A4_data = {
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
    };
