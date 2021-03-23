function [r1WithSchroedingerEquation] = ...
    calculateR1WithSchroedingerEquation(theta,phi ...
    ,nearestNeighbourDistancesPow3,deltaT,path2ConstantsFile)


%% Define Constants
constants = readConstantsFile(path2ConstantsFile);

Nm = constants.Nm;
hbar = constants.hbar;
gammaRad = constants.gammaRad;
B0 = constants.B0;
mu0 = constants.mu0;
w0 = gammaRad*B0;
alpha = mu0/(4*pi)*gammaRad^2*hbar^2/Nm^3;

%% Define Simulation Paramter
[numberOfHs,numberOfTimeSteps] = size(theta);
timeAxis = linspace(deltaT,numberOfTimeSteps*deltaT,numberOfTimeSteps);

% Pauli matrices
pauliSpinX = 1/2*[[0 1]' [1 0]'];
pauliSpinY = 1/2*[[0 -1i]' [1i 0]'];
pauliSpinZ = 1/2*[[1 0]' [0 -1]'];

pauliSpinCreation = (pauliSpinX+1i*pauliSpinY);
pauliSpinAnnihilation = (pauliSpinX-1i*pauliSpinY);

% Spin Matrices for a two spin System
identityMatrix = eye(2,2);

firstSpinZ=kron(pauliSpinZ,identityMatrix);
firstSpinCreation=kron(pauliSpinCreation,identityMatrix);
firstSpinAnnihilation=kron(pauliSpinAnnihilation,identityMatrix);

secondSpinZ=kron(identityMatrix,pauliSpinZ);
secondSpinCreation=kron(identityMatrix,pauliSpinCreation);
secondSpinAnnihilation=kron(identityMatrix,pauliSpinAnnihilation);

% Operators from interaction
A = secondSpinZ*firstSpinZ;
B = -1/4*(secondSpinCreation*firstSpinAnnihilation ...
    +secondSpinAnnihilation*firstSpinCreation);
C = -(3/2)*(secondSpinCreation*firstSpinZ ...
    +secondSpinZ*firstSpinCreation);
D = -(3/2)*(secondSpinAnnihilation*firstSpinZ ...
    +secondSpinZ*firstSpinAnnihilation);
E = -(3/4)*(secondSpinCreation*firstSpinCreation);
F = -(3/4)*(secondSpinAnnihilation*firstSpinAnnihilation);

% Spherical Harmonics
SH_A = (1-3*cos(theta).^2)./nearestNeighbourDistancesPow3;
SH_B = (1-3*cos(theta).^2)./nearestNeighbourDistancesPow3;
SH_C = sin(theta).*cos(theta).*exp(-1i*phi)./nearestNeighbourDistancesPow3;
SH_D = sin(theta).*cos(theta).*exp(1i*phi)./nearestNeighbourDistancesPow3;
SH_E = sin(theta).^2.*exp(-2i*phi)./nearestNeighbourDistancesPow3;
SH_F = sin(theta).^2.*exp(2i*phi)./nearestNeighbourDistancesPow3;

% Possible states
psi_pp=kron([1 0],[1 0])';
psi_mp=kron([0 1],[1 0])';
psi_pm=kron([1 0],[0 1])';
psi_mm=kron([0 1],[0 1])';

% Transition rates
transitionRates_mm2pm=zeros(numberOfTimeSteps,1);
transitionRates_mm2mp=zeros(numberOfTimeSteps,1);
transitionRates_mm2pp=zeros(numberOfTimeSteps,1);

% Transition probabilities
transitionProbabilities_mm2pm=zeros(numberOfTimeSteps,1);
transitionProbabilities_mm2mp=zeros(numberOfTimeSteps,1);
transitionProbabilities_mm2pp=zeros(numberOfTimeSteps,1);

%% Simulation of interacting particles
r1 = zeros(1,numberOfHs);
for atomNumber = 1:numberOfHs
    
    disp(['Solve Schroedinger Equation for Atom Number ' ...
        ,num2str(atomNumber),' of ',num2str(numberOfHs)])
    H0 = -w0*hbar*(firstSpinZ+secondSpinZ);
    psi = psi_mm;
    
    for timeStep = 1:numberOfTimeSteps
        
        H1 = alpha*(A*SH_A(atomNumber,timeStep) ...
            +B*SH_B(atomNumber,timeStep) ...
            +C*SH_C(atomNumber,timeStep)...
            +D*SH_D(atomNumber,timeStep) ...
            +E*SH_E(atomNumber,timeStep) ...
            +F*SH_F(atomNumber,timeStep));
        
        U0 = expm(-1i*H0*timeAxis(timeStep)/hbar);
        H_IP = U0'*H1*U0;
        
        % Euler
        psi = psi+H_IP*psi*deltaT/(1i*hbar);
        psi = psi/norm(psi);
        
        % calculate transition probabilities
        transitionProbabilities_mm2pm(timeStep)=abs(dot(psi,psi_pm))^2;
        transitionProbabilities_mm2mp(timeStep)=abs(dot(psi,psi_mp))^2;
        transitionProbabilities_mm2pp(timeStep)=abs(dot(psi,psi_pp))^2;
        
        % calculate transition rates
        transitionRates_mm2pm(timeStep) = abs(dot(psi,psi_pm))^2 ...
            /(deltaT*timeStep);
        transitionRates_mm2mp(timeStep) = abs(dot(psi,psi_mp))^2 ...
            /(deltaT*timeStep);
        transitionRates_mm2pp(timeStep) = abs(dot(psi,psi_pp))^2 ...
            /(deltaT*timeStep);
    end
    
    %R1 = 2*W_(singleFlip)+2*W_(doubleFlip)
    r1(atomNumber) = transitionRates_mm2pm(end) ...
        +transitionRates_mm2mp(end)+2*transitionRates_mm2pp(end);
end

r1WithSchroedingerEquation = sum(r1);

end

%% Runge Kutta-------
%         k1=SG(Psi,H_IP);
%         k2=SG(Psi+k1*DeltaT/2,H_IP);
%         k3=SG(Psi+k2*DeltaT/2,H_IP);
%         k4=SG(Psi+k3*DeltaT,H_IP);
%         Psi=Psi+1/6*(k1+2*k2+2*k3+k4)*DeltaT;
%         Psi=Psi/norm(Psi);
%%---------------------------


