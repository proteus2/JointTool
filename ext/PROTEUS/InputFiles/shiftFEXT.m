function fext = shiftFEXT(fext,xyzshift,str)

global GRAPHFLAG

% --- Shift external masses based on shift of wing data
for i=1:length(fext.type)
    fext.location{i}(:,1) = fext.location{i}(:,1)-xyzshift(1);
    fext.location{i}(:,2) = fext.location{i}(:,2)-xyzshift(2);
    fext.location{i}(:,3) = fext.location{i}(:,3)-xyzshift(3);
end
% ---

% --- Process lumped data
fext = fext_inp(str,fext);
% ---

% ---
if GRAPHFLAG == 1
    figure(1); hold on
    for ifext = 1 : length(fext.type)
        for jfext = 1 : size(fext.magnitude{ifext},1)
           plot3(fext.location{ifext}(jfext,1),fext.location{ifext}(jfext,2),fext.location{ifext}(jfext,3),'blacks','MarkerFaceColor','none','LineStyle','-')
        end
    end
end