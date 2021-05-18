clc
clear all

  clc
  clear all

path2Results= ['C:\Users\maxoe\Google Drive\Promotion\BackUps\' ...
    'CodeBackUpFelix_OrientDep\Results\'];
load(['C:\Users\maxoe\Google Drive\Promotion\Data\Myelin\' ...
    'Water_H_50ns_025ps_nH100.mat']);
Shift=1;
water_H1=Water_H_50ns_025ps_nH100(1:2:end,:,10000:Shift:end);
water_H2=Water_H_50ns_025ps_nH100(2:2:end,:,10000:Shift:end);

if size(water_H1,1) ~= size(water_H2)
    water_H1 = water_H1(1:end-1,:,:);
end

X=squeeze(water_H1(:,1,1:end))-squeeze(water_H2(:,1,1:end));
Y=squeeze(water_H1(:,2,1:end))-squeeze(water_H2(:,2,1:end));
Z=squeeze(water_H1(:,3,1:end))-squeeze(water_H2(:,3,1:end));

distances = sqrt(X.^2+Y.^2+Z.^2);

for testCase = 1:20
   randomAtomNumber = randi(size(water_H1,1));
   disp('=========================')
   disp(['X = ' num2str(X(randomAtomNumber,100000))])
   disp(['Distances = ' num2str(distances(randomAtomNumber,1))])
   disp(['H1 = ' num2str(water_H1(randomAtomNumber,100000))])
   
   
end


figure(1)
moleculePositions = Water_H_50ns_025ps_nH100(:,:,100000);
plot3(moleculePositions(:,1),moleculePositions(:,2),moleculePositions(:,3),'*')

figure(2)
atomPositionsFirstH = water_H1(:,:,100000);
atomPositionsSecondH = water_H2(:,:,100000);
plot3(atomPositionsFirstH(:,1),atomPositionsFirstH(:,2),atomPositionsFirstH(:,3),'+')
hold on
plot3(atomPositionsSecondH(:,1),atomPositionsSecondH(:,2),atomPositionsSecondH(:,3),'*')
plot3(X(:,100000),Y(:,100000),Z(:,100000),'o')
hold off
grid on
legend('first H', 'second H','Position Calculation')






