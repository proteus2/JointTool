function lumped = shiftNSM(lumped,xyzshift,str)

global GRAPHFLAG

% if isfield(lumped,'type')
%     LastEntryLumped = length(lumped.type);
% else
    LastEntryLumped = 1;
% end

% --- Shift external masses based on shift of wing data
for i=1:length(lumped.location)-LastEntryLumped
    lumped.location{LastEntryLumped+i}(:,1) = lumped.location{LastEntryLumped+i}(:,1)-xyzshift(1);
    lumped.location{LastEntryLumped+i}(:,2) = lumped.location{LastEntryLumped+i}(:,2)-xyzshift(2);
    lumped.location{LastEntryLumped+i}(:,3) = lumped.location{LastEntryLumped+i}(:,3)-xyzshift(3);
end
% ---

% --- Process lumped data
lumped = mext_inp(str,lumped);
% ---

% ---
if GRAPHFLAG == 1
    figure(1); hold on
%     for iNSmass = 1 : length(lumped.location)-LastEntryLumped
%            plot3(lumped.location{LastEntryLumped+iNSmass}(:,1),lumped.location{LastEntryLumped+iNSmass}(:,2),lumped.location{LastEntryLumped+iNSmass}(:,3),'redo','MarkerFaceColor','red','MarkerSize',6)
%     end
    for iNSmass = 1 : length(lumped.location)
           plot3(lumped.location{iNSmass}(:,1),lumped.location{iNSmass}(:,2),lumped.location{iNSmass}(:,3),'redo','MarkerFaceColor','red','MarkerSize',6)
    end
end