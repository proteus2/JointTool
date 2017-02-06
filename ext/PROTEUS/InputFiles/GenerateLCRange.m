function [lcrange] = GenerateLCRange(xlsFileName,model)

LCRangeExcel = xlsread(xlsFileName,'Loadcase Selection');

lcase.nz=LCRangeExcel(1,1);
lcase.trim=LCRangeExcel(2,1);
lcase.alpha0=LCRangeExcel(3,1);
lcase.gustflag=LCRangeExcel(4,1);
lcase.v=(LCRangeExcel(9,1):LCRangeExcel(10,1):LCRangeExcel(11,1))';
lcase.alt=(LCRangeExcel(12,1):LCRangeExcel(13,1):LCRangeExcel(14,1))';
if lcase.gustflag == 1
    if ~isnan(LCRangeExcel(16,1))
        lcase.gustlength=(LCRangeExcel(16,1):LCRangeExcel(17,1):LCRangeExcel(18,1))';
    else
        lcase.gustlength = [];
    end
    if ~isnan(LCRangeExcel(16,2))
        lcase.gustlength=[lcase.gustlength;(LCRangeExcel(16,2):LCRangeExcel(17,2):LCRangeExcel(18,2))'];
    end
end
   
EnvExcel = xlsread(xlsFileName,'Flight Envelope');

envelope.velocity = EnvExcel(:,1);
envelope.altitude = EnvExcel(:,2);

lcenvelope=[];
for i=1:length(lcase.v)
    for j=1:length(lcase.alt)
        [in,edge] = inpolygon(lcase.v(i),lcase.alt(j),envelope.velocity,envelope.altitude);
        if in==1 || edge==1
            [rho,a] = stdatmo(lcase.alt(j));     
            lcenvelope(end+1,1:5)=[lcase.v(i),lcase.alt(j),0.5*stdatmo(lcase.alt(j))*lcase.v(i)^2,tas2eas(lcase.v(i),lcase.alt(j)),lcase.v(i)/a;];
        end                    
    end
end                
            
for i=1:size(LCRangeExcel,2)
    fuel_level(:,i)=(LCRangeExcel(6,i):LCRangeExcel(7,i):LCRangeExcel(8,i))';
end

lcmat=repmat(lcenvelope,size(fuel_level,1),1);

fuel_level = fuel_level(reshape(repmat((1:size(fuel_level,1))',1,size(lcenvelope,1))',[],1),:);
            
lcrange.fuel_level=mat2cell(fuel_level,ones(size(fuel_level,1),1));
nlc = size(lcmat,1);
lcrange.nz = lcase.nz*ones(nlc,1);
lcrange.gustflag = lcase.gustflag*ones(nlc,1);
lcrange.trim = lcase.trim*ones(nlc,1);
lcrange.alpha0 = lcase.alpha0*ones(nlc,1);
lcrange.EAS = lcmat(:,4);
lcrange.M = lcmat(:,5);
lcrange.H = lcmat(:,2);
lcrange.LCload = NaN*ones(nlc,1);

lcrange.gustlength = lcase.gustlength;
