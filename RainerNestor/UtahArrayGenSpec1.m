% parameters
NumElecX = 10;
NumElecY = 10;

ElecLength = 1.5;
ElecDistance = 0.4;
BasePlateThickness = 0.5;

outputfilepath = './Utah10by10.ecp';

% calculated
NumElec = NumElecX * NumElecY;

ArrayLengthX = ElecDistance * (NumElecX-1);
ArrayLengthY = ElecDistance * (NumElecY-1);

Array.x = zeros(NumElecX,NumElecY);
Array.y = zeros(NumElecX,NumElecY);
Array.z = ones(NumElecX,NumElecY) * (- ElecLength);

for thisY = 1:NumElecY
    Array.x(:,thisY) = linspace(0,ArrayLengthX,NumElecX);
end

for thisX = 1:NumElecX
    Array.y(thisX,:) = linspace(0,ArrayLengthY,NumElecY);
end

Array.x = Array.x - ArrayLengthX/2;
Array.y = Array.y - ArrayLengthY/2;

% generate .ecp file
fp = fopen(outputfilepath,'w+');
fprintf(fp,'FileVersion: 1\n');
fprintf(fp,'AlongElectrodesAxis: 2\n');
fprintf(fp,'Origin:              0.0 0.0 0.0\n');
fprintf(fp,'NrOfContactPoints:   %d\n',NumElec);

for thisX = 1:NumElecX
    for thisY = 1:NumElecY
        fprintf(fp,'%.2f %.2f %.2f\n',Array.x(thisX,thisY),Array.y(thisX,thisY),Array.z(thisX,thisY));
    end
end
fclose(fp);