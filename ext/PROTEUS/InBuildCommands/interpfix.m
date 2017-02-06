function [datay] = interpfix(datax,N)

elml = (datax(end)-datax(1))/(N-1);
datay=[];
for i=1:length(datax)-1
    numelloc = ceil((datax(i+1)-datax(i))/elml);
    datay = [datay;linspace(datax(i),datax(i+1),numelloc+1)'];
end
datay = unique(datay);