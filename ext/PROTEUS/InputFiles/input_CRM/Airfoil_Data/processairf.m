clear all
close all
clc

for i=1:16
   airf_raw{i} = dlmread(['airfoil',num2str(i),'.txt']); 
    
   [~,leind] = min(airf_raw{i}(:,1));
   
   leloc = airf_raw{i}(leind,:);
   teloc = (airf_raw{i}(1,:)+airf_raw{i}(end,:))/2;
   
   twist(i,1) = -atand((teloc(3)-leloc(3))/(teloc(1)-leloc(1)));
   
   Rrot = [cosd(twist(i,1)) ,0 ,-sind(twist(i,1));
       0,1,0;
       sind(twist(i,1)),0,cosd(twist(i,1))];
   for j=1:size(airf_raw{i},1)
        airf_rot{i}(j,:) = (Rrot*(airf_raw{i}(j,:)-leloc)')';
        
        
   end
   airf_fin_top{i}(:,1) = airf_rot{i}(leind:-1:1,1)/max(airf_rot{i}(:,1));
   airf_fin_top{i}(:,2) = airf_rot{i}(leind:-1:1,3)/max(airf_rot{i}(:,1));
   airf_fin_bot{i}(:,1) = airf_rot{i}(leind:end,1)/max(airf_rot{i}(:,1));
   airf_fin_bot{i}(:,2) = airf_rot{i}(leind:end,3)/max(airf_rot{i}(:,1));
   yloc(i,1) = 0.0254*mean(airf_raw{i}(:,2));
   
   dlmwrite(['airfoil',num2str(i),'_scaled.txt'],[airf_fin_top{i},airf_fin_bot{i}],'delimiter',' ','precision','%6f');
end


figure
hold all
airf = 7;
% plot(airf_raw{14}(:,1),airf_raw{14}(:,3))
% plot(airf_rot{14}(:,1),airf_rot{14}(:,3))
plot(airf_fin_top{airf}(:,1),airf_fin_top{airf}(:,2))
plot(airf_fin_bot{airf}(:,1),airf_fin_bot{airf}(:,2))
plot(airf_fin_bot{airf}(:,1),(airf_fin_top{airf}(:,2)+airf_fin_bot{airf}(:,2))/2)
% axis equal