function calculateR1RateRAMandPERFORMANCEoptimized(configuration)

[path2Data,path2Save,path2ConstantsFile] = setUpSystemBasedOnMachine( ...
    configuration);

[trajectoryX,trajectoryY,trajectoryZ] = loadTrajectoriesFromData( ...
    configuration,path2Data);
fileName = configuration.fileName; %#ok<NASGU>

clearvars -except  trajectoryX trajectoryY trajectoryZ configuration ...
    fileName path2Save path2ConstantsFile


%% Define constants
disp('Defining constants')
deltaT = configuration.deltaT;
constants = readConstantsFile(path2ConstantsFile);

picoSecond = constants.picoSecond; %#ok<NASGU>
hbar = constants.hbar;
gammaRad = constants.gammaRad;
B0 = configuration.B0;
mu0 = constants.mu0;
Nm = constants.Nm;
DD = 3/4*(mu0/(4*pi)*hbar*gammaRad^2)^2/(Nm^6);
omega0 = gammaRad*B0;

%% Define simulation parameters
disp('Defining simulation parameters')
[numberOfHs,timeSteps] = size(trajectoryX);
lags = round(configuration.fractionForLags*timeSteps);
nearestNeighbours = configuration.nearestNeighbours;
atomsToCalculate = configuration.atomsToCalculate;
startDateOfSimulation = datestr(now,'yyyymmdd'); %#ok<NASGU>

%% Start simulation
stopWatch = zeros(1,atomsToCalculate);

orientationAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.fibreOrientations,';','numeric'));
positionAngles = deg2rad(getValuesFromStringEnumeration( ...
    configuration.myelinPositions,';','numeric'));

fibreOrientationsCount = size(orientationAngles,2);
positionsInMyelinCount = size(positionAngles,2);

meanPositions = single([mean(trajectoryX,2) mean(trajectoryY,2) ...
    mean(trajectoryZ,2)]); %#ok<NASGU>

atomCounter = 0;
atomIndex = zeros(1,atomsToCalculate);

r1WithPerturbationTheory = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate);

nearestNeighboursX = zeros(nearestNeighbours,size(trajectoryX,2));
nearestNeighboursY = zeros(nearestNeighbours,size(trajectoryX,2));
nearestNeighboursZ = zeros(nearestNeighbours,size(trajectoryX,2));
nearestNeighbourDistancesPow3 = zeros(size(nearestNeighboursX));

distances = zeros(size(trajectoryX));
inverseMeanDistances = 1./mean(distances,2);

rotatedX = zeros(size(nearestNeighboursX));
rotatedY = zeros(size(nearestNeighboursX));
rotatedZ = zeros(size(nearestNeighboursX));


polarAngle = zeros(size(rotatedX));
azimuthAngle = zeros(size(rotatedX));
hypotuseXY = zeros(size(rotatedX));

firstOrderSphericalHarmonic = zeros(size(rotatedX,1),size(rotatedY,2));
secondOrderSphericalHarmonic = zeros(size(rotatedX,1),size(rotatedY,2));
[nearestNeighbours,timeSteps] = size(firstOrderSphericalHarmonic);
zeroPaddingLength = 2^(nextpow2(timeSteps)+1);

correlationFunction1W0 = zeros(nearestNeighbours,zeroPaddingLength);
correlationFunction2W0 = zeros(nearestNeighbours,zeroPaddingLength);

shiftForCorrelationFunction = configuration.shiftForCorrelationFunction;
tmp = false(1,lags);
tmp = tmp(1:shiftForCorrelationFunction:end);
correlationFunction1W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,length(tmp) ...
    ,'like',single(1j));
correlationFunction2W0Saver = zeros(fibreOrientationsCount ...
    ,positionsInMyelinCount,atomsToCalculate,length(tmp) ...
    ,'like',single(1j));
clearvars tmp



disp('Starting Simulation.')

for atomNumber = randperm(numberOfHs)
    wholeTic = tic;
    transformationTic = tic;
    disp('=============================')
    disp(['Atom number: ' num2str(atomNumber)])
    disp(['Calculated ' num2str(atomCounter) ' atoms.'])
    atomCounter = atomCounter+1;
    atomIndex(atomCounter) = atomNumber;
    
    calculateRelativePositions();
    
    findNearestNeighbours();
    
    clear relativeX relativeY relativeZ
    
    transformationTime = toc(transformationTic);
    disp(['Time needed for transforming to relative variables: ' ...
        num2str(transformationTime) 's']);
    positionsTic = tic;
    for positionNumber = 1:positionsInMyelinCount
        positionAngle = positionAngles(positionNumber);
        zAxis = [0 0 1];
        rotationMatrixPosition = get3DRotationMatrix( ...
            positionAngle,zAxis);
        orientationsTic = tic;
        for orientationNumber = 1:fibreOrientationsCount
            orientationAngle = orientationAngles(orientationNumber);
            yAxis = [0 1 0];
            rotationMatrixOrientation = get3DRotationMatrix( ...
                orientationAngle,yAxis);
            rotationMatrix = ...
                rotationMatrixOrientation*rotationMatrixPosition;
            
            rotateTrajectoriesWithRotationMatrix();
            
            calculateR1RatesAndCorrelationFunctions();
            
        end
        disp(['Time needed for all orientations at position ' ...
            num2str(rad2deg(positionAngle)) ': ' ...
            num2str(toc(orientationsTic)) 's'])
    end
    disp(['Time needed for all variations : ' ...
        num2str(toc(positionsTic)) 's'])
    if mod(atomCounter,5) == 0
        disp('Start saving data.')
        savingTic = tic;
        lastSavingDate = datestr(now,'yyyymmdd_HHMM'); %#ok<NASGU>
        save(path2Save ,'r1WithPerturbationTheory' ...
            ,'correlationFunction1W0Saver' ...
            ,'correlationFunction2W0Saver' ...
            ,'meanPositions','deltaT' ...
            ,'timeSteps','lags','atomCounter','B0' ...
            ,'orientationAngles','positionAngles' ...
            ,'atomIndex','nearestNeighbours','atomsToCalculate' ...
            ,'fileName','stopWatch','numberOfHs' ...
            ,'startDateOfSimulation','lastSavingDate','constants' ...
            ,'shiftForCorrelationFunction','DD','-v7.3')
        disp(['Time needed for saving variables: ' ...
            num2str(toc(savingTic)) 's'])
    end
    
    stopWatch(atomCounter) = toc(wholeTic);
    disp(['Atom ' num2str(atomNumber) ' done.'])
    disp(['Overall time needed for one Atom: ' ...
        num2str(stopWatch(atomCounter)) 's'])
end

    function calculateRelativePositions()
        
        trajectoryX = bsxfun(@minus,trajectoryX,trajectoryX(atomNumber,:));
        trajectoryY = bsxfun(@minus,trajectoryY,trajectoryY(atomNumber,:));
        trajectoryZ = bsxfun(@minus,trajectoryZ,trajectoryZ(atomNumber,:));
        
    end

    function findNearestNeighbours()
        
        distances = sqrt(trajectoryX.^2+trajectoryY.^2+trajectoryZ.^2);
        inverseMeanDistances = 1./mean(distances,2);
        inverseDistances = 0;
        sumIndex = 1;
        
        for i=1:nearestNeighbours+1
            smallestDistance = max(inverseMeanDistances);
            closestNeighbourId = inverseMeanDistances == smallestDistance;
            numberOfClosestElements = sum(closestNeighbourId);
            inverseDistances(sumIndex:sumIndex+numberOfClosestElements-1) = ...
                smallestDistance;
            nearestNeighbourIndex(sumIndex:sumIndex+numberOfClosestElements-1) = ...
                find(inverseMeanDistances == smallestDistance);
            sumIndex = sumIndex + numberOfClosestElements;
            inverseMeanDistances(closestNeighbourId) = min(inverseMeanDistances);
        end
        
        neighbourIds = false(size(distances));
        neighbourIds(nearestNeighbourIndex) = true;
        neighbourIds(atomNumber) = false;
        
        nearestNeighboursX = trajectoryX(neighbourIds,:);
        nearestNeighboursY = trajectoryY(neighbourIds,:);
        nearestNeighboursZ = trajectoryZ(neighbourIds,:);
        nearestNeighbourDistancesPow3 = distances(neighbourIds,:).^3;
    end

    function rotateTrajectoriesWithRotationMatrix()
        
        rotatedX = (rotationMatrix(1,1)*nearestNeighboursX ...
            +rotationMatrix(1,2)* nearestNeighboursY ...
            +rotationMatrix(1,3)*nearestNeighboursZ );
        rotatedY = (rotationMatrix(2,1)*nearestNeighboursX ...
            +rotationMatrix(2,2)* nearestNeighboursY ...
            +rotationMatrix(2,3)*nearestNeighboursZ );
        rotatedZ = (rotationMatrix(3,1)*nearestNeighboursX ...
            +rotationMatrix(3,2)* nearestNeighboursY ...
            +rotationMatrix(3,3)*nearestNeighboursZ );
        
    end

    function calculateR1RatesAndCorrelationFunctions()
        
        spectralDensity1W0 = single(0);
        spectralDensity2W0 = single(0);
        
        
        transformToSphericalCoordinates();
        function transformToSphericalCoordinates()
            
            hypotuseXY = hypot(rotatedX,rotatedY);
            polarAngle = pi/2-atan2(rotatedZ,hypotuseXY);
            azimuthAngle = atan2(rotatedY,rotatedX);
            
        end
        
        calculateSphericalHarmonics();
        function calculateSphericalHarmonics()
            
            firstOrderSphericalHarmonic = sin(polarAngle).*cos(polarAngle) ...
                .*exp(-1i*azimuthAngle)./nearestNeighbourDistancesPow3;
            secondOrderSphericalHarmonic = sin(polarAngle).^2 ...
                .*exp(-2i*azimuthAngle)./nearestNeighbourDistancesPow3;
        end
        
        calculateCrossCorrelationFunction();
        function calculateCrossCorrelationFunction()
            
            firstOrderSphericalHarmonic = fft(firstOrderSphericalHarmonic ...
                ,zeroPaddingLength,2);
            % Convolution
            correlationFunction1W0 = ifft(firstOrderSphericalHarmonic ...
                .*conj(firstOrderSphericalHarmonic),[],2)/timeSteps;
            correlationFunction1W0 = sum(correlationFunction1W0(:,1:lags));
            
            secondOrderSphericalHarmonic = fft(secondOrderSphericalHarmonic ...
                ,zeroPaddingLength,2);
            % Convolution
            correlationFunction2W0 = ifft(secondOrderSphericalHarmonic ...
                .*conj(secondOrderSphericalHarmonic),[],2)/timeSteps;
            correlationFunction2W0 = sum(correlationFunction2W0(:,1:lags));
        end
        
        writeCorrelationFunctionSaver();
        function writeCorrelationFunctionSaver()
            correlationFunction1W0Saver(orientationNumber ...
                ,positionNumber,atomCounter,:) = single( ...
                correlationFunction1W0( ...
                1:shiftForCorrelationFunction:end)); %#ok<SETNU>
            correlationFunction2W0Saver(orientationNumber ...
                ,positionNumber,atomCounter,:) = single( ...
                correlationFunction2W0( ...
                1:shiftForCorrelationFunction:end)); %#ok<SETNU>
        end 
        
        calculateSpectralDensities();
        function calculateSpectralDensities()
            
            spectralDensity1W0 = 2*(deltaT*cumsum(correlationFunction1W0 ...
                .*exp(-1i*omega0*deltaT*(0:lags-1))));
            spectralDensity2W0 = 2*(deltaT*cumsum(correlationFunction2W0 ...
                .*exp(-2i*omega0*deltaT*(0:lags-1))));
            
            spectralDensity1W0 = spectralDensity1W0(end);
            spectralDensity2W0 = spectralDensity2W0(end);
            
        end
        
        calculateR1WithSpectralDensity();
        function calculateR1WithSpectralDensity()
            
            r1WithPerturbationTheory(orientationNumber,positionNumber ...
                ,atomCounter) = DD*3/2*( ...
                abs(real(mean(spectralDensity1W0))) ...
                +abs(real(mean(spectralDensity2W0)))); %#ok<SETNU>
            
        end
        
    end

end


