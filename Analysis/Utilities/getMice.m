function category = getMice(fileName)

numPart = regexp(fileName, '\d{6}_\d{2}_\d{2}_\d{2}', 'match', 'once');

Jeremy_names = {
    '191018_13_39_41';
    '191018_13_56_55';
    '191018_14_30_00';
    '191018_14_11_33';
    };
Bernie_names = {
   '191209_13_44_12';
   '191209_14_04_14';
   '191209_14_32_39';
   '191209_15_01_22';
   '191209_14_18_13';
   '191209_14_46_58';
    };
Bill_names  = {
    '200130_13_21_13';
    '200130_13_36_14';
    '200130_13_49_09';
    '200130_14_02_12';
    '200130_14_15_24';
    '200130_14_29_30';
    };
Nigel_names = {
    '171212_16_19_37';
    };

if ismember(numPart, Jeremy_names)
    category = 'Jeremy';
elseif ismember(numPart, Bernie_names)
    category = 'Bernie';
elseif ismember(numPart, Bill_names)
    category = 'Bill';
elseif ismember(numPart, Nigel_names)
    category = 'Nigel';
else
    category = 'Unknown';
end

end
