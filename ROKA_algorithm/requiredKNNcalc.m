function [reqKnn] = requiredKNNcalc(PC,ptCloud,userSPr)
rn=randi([1 length(PC(:,1))],1,10000);
for i = 1:length(rn)
[idxKrn, ptdistrn] = findNearestNeighbors(ptCloud,PC(rn(i),1:3),5);
meandrn(i)=mean(ptdistrn(2:end));
end
Arn=pi *mean(meandrn)^2;
ptRho=5/Arn;

Auser=pi * userSPr^2;
reqKnn=ceil(Auser*1.65*(ptRho));
end
