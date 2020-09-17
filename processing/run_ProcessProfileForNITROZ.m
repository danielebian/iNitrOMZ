% 
addpath('/Users/yangsi/Box Sync/UCLA/Data/1-Dmodel/Colette');
File1Path = '/Users/yangsi/Box Sync/UCLA/Data/1-Dmodel/CompilationOct232018/190606_offshore_station.xlsx';
File2Path = '/Users/yangsi/Box Sync/UCLA/Data/1-Dmodel/CompilationOct232018/190606_coastal_station.xlsx';

ColumnNames={'Depth','Pressure','Longitude','Latitude','temperature','salinity','sigmat','o2','no3','no2','po4','nh4','d15no3','d15no2','n2o','d15Na','d15Nb','d15NBulk'};
compilation_offshore =  ProcessProfileForNITROZ(File1Path, ColumnNames);
compilation_coastal =  ProcessProfileForNITROZ(File2Path, ColumnNames);

save('/Users/yangsi/Box Sync/UCLA/Data/1-Dmodel/Colette/compilation_offshore.mat','compilation_offshore');
save('/Users/yangsi/Box Sync/UCLA/Data/1-Dmodel/Colette/compilation_coastal.mat','compilation_coastal');
