function [A] = globel(Al,eft)

A=zeros(max(max(eft)));
for i=1:size(eft,1)
   A(eft(i,:),eft(i,:))= A(eft(i,:),eft(i,:))+squeeze(Al(i,:,:));
end