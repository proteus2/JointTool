function loadcase=GenerateROMLoadcase(ExcelFile)

LoadcaseExcel = xlsread(ExcelFile,'ROM Loadcase');

loadcase.M  (1,1)  = LoadcaseExcel(1,2);
loadcase.EAS(1,1)  = LoadcaseExcel(1,3);
loadcase.H  (1,1)  = LoadcaseExcel(1,4);
loadcase.LCload(1,1) = LoadcaseExcel(1,5);
loadcase.nz (1,1)  = LoadcaseExcel(1,6);
loadcase.gustflag(1,1) = LoadcaseExcel(1,7);
loadcase.trim(1,1) = LoadcaseExcel(1,8);
loadcase.alpha0(1,1) = LoadcaseExcel(1,9);
loadcase.fuel_level{1} = LoadcaseExcel(1,10:end);

