function [Xall,Tall,Lall,A_nodes,Ma,Mb] = DR(Assembly)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Define Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Links = Assembly.Elements2Dof;
X0 = Assembly.Positions;
Type = Assembly.ElementType;
Supports = Assembly.Supports;
A = Assembly.ElementArea;
Young = Assembly.ElementYoung;
L0 = Assembly.ElementInitialLength;
A_a0 = Assembly.ElementInitialAnglea;
A_b0 = Assembly.ElementInitialAngleb;
P = Assembly.Prestress;
dL = Assembly.ElementLengthChange;
Iy = Assembly.ElementInertia;

nCase = size(Assembly.ExternalLoad,2);
nElements = size(Links,1);		% number of links
nDim = size(Links,2)/2;                  % problem dimension
nNodes = size(X0,1)/nDim;			% node number
DeltaT = 0.01;                     % time interval (small value)

A_nodes = zeros(nNodes,1);          %Initial orientation of the nodes

iSupports = ~Supports;

% INITIALIZE OUTPUT
Xall = zeros(nNodes*nDim,nCase); Tall = zeros(nElements,nCase); Lall = zeros(nElements,nCase);

% FOR ALL LOAD CASES
for i = 1:nCase
    
X = X0;                         % initial coordinates
    
% EXTERNAL CHANGES
Fe = Assembly.ExternalLoad(:,i);
Lnew = L0 + dL(:,i);

% COMPUTING OF RESIDUAL FORCES
[RF,RM,L,T,Ma,Mb] = DrResiduals(Fe,Links,Type,X,Young,A,Iy,Lnew,A_nodes,A_a0,A_b0,P);

% Calling the function for the computing of masses
M = DrMasses(Links,Type,X,Supports,Young,A,L,T,DeltaT);
I = zeros(size(M,1)/2,1);
iSupportsM = zeros(size(M,1)/2,1);
for i=1:size(M,1)/2
    
        I(i,1)=M(2*i-1)*1000;                          %Augmenter évantuelement
        iSupportsM(i,1)=iSupports(2*i-1);
    
end

% COMPUTING OF INITIAL VELOCITIES
V = 0.5*DeltaT*rdivide(RF,M).*iSupports;
OM = 0.5*DeltaT*rdivide(RM,I).*iSupportsM;

% DYNAMIC RELAXATION: ALGORITHM IMPLEMENTATION
KE = 0;                       % kinetic energy set to zero
nbCycle = 0;					% number of cycles
nbReset = 0;					% number of resets

% main loop of the dynamic relaxation algorithm
while 1	% while not converged
    
    % velocity computation I
    V = (V + (DeltaT*rdivide(RF,M))).*iSupports;
    X = X + DeltaT*V;
    Y = norm(DeltaT*V);
    
    % angular velocity computation
    OM = (OM + (DeltaT*rdivide(RM,I))).*iSupportsM;
    A_nodes = A_nodes + DeltaT*OM;
    
    % kinetic energy calculation
    KEdt = sum(0.5*(V.*V).*M);
    
    % kinetic damping computation
    if KEdt<=KE          % KE Peak detected : if (current kinetic energy) < (previous kinetic)
        X = X - (1.5*DeltaT*V) + (0.5*DeltaT*DeltaT*(RF./M));
        [RF,RM,~] = DrResiduals(Fe,Links,Type,X,Young,A,Iy,Lnew,A_nodes,A_a0,A_b0,P);
        % velocity computation II
        V = 0.5*DeltaT*(RF./M).*iSupports;
        OM = 0.5*DeltaT*(RM./I).*iSupportsM;
        
        RF = RF.*iSupports;
        KEdt = 0;
        nbReset = nbReset + 1 ;
        
        % Convergence verification through calculating Residual Norm
        NormRes = norm(RF);
        if NormRes < 0.0001 || nbCycle > 5000 || Y < 0.0001
            break
        end
    end
    
    % Recalculate Residuals, Masses and update Current kinetic energy
    [RF,RM,L,T,Ma,Mb] = DrResiduals(Fe,Links,Type,X,Young,A,Iy,Lnew,A_nodes,A_a0,A_b0,P);    % calling the function for the computing of residual forces
    [M] = DrMasses(Links,Type,X,Supports,Young,A,L,T,DeltaT);               % calling the function for the computing of masses
    
    I = zeros(size(M,1)/2,1);
    iSupportsM = zeros(size(M,1)/2,1);
    for i=1:size(M,1)/2
    
        I(i,1)=M(2*i-1)*10;                          %Augmenter évantuelement
        iSupportsM(i,1)=iSupports(2*i-1);
    
    end
    
    KE = KEdt; % current kinetic energy is stored as previous one
    
   
    nbCycle = nbCycle + 1;
    
    plot(A_nodes)

    
%     For Dynamic animation    
%     Xall(:,nbCycle) = X;
%     Tall(:,nbCycle) = T;
%     Lall(:,nbCycle) = L;
    
end

Xall(:,i) = X;
Tall(:,i) = T;
Lall(:,i) = L;

end

function [M] = DrMasses(Links,Type,X,Bound,Young,A,L,T,DeltaT)

% Reset Residuals to nodal forces values
nElements = size(Links,1);					% number of links
nDim = size(Links,2)/2;                  % problem dimension
BoundMass = 1.e100;      				% large value for boundary nodes
minM = 0.005;							% small value for avoiding zero-masses
amplMass = 1;							% amplification of fictitious nodal masses
DT2 = DeltaT*DeltaT;						% square of the time interval

% Mass Reset
M = X*0 + Bound*BoundMass;

EA = Young.*A;
L2 = L.^2;

k = (EA+T.*Type)./L; % axial rigidity 

% nodal mass computation
for i=1:nElements							% for all links
    for j = 1:nDim
        DX = abs(X(Links(i,j+nDim))-X(Links(i,j)));		% difference in coordinates
        kX = k(i)*((DX*DX)/L2(i)); % axial rigidity for every direction		difference in coordinates ^2 / length ^2 = cosinus of direction

        MX = kX*2*DT2*amplMass + minM;		% mass calculation                          mass in direction = (axial rigidity on the direction * TimeInterval^2 * amplification) + min mass
        M(Links(i,j)) = M(Links(i,j)) + MX;				% mass for node A
        M(Links(i,j+nDim)) = M(Links(i,j+nDim)) + MX;	% mass for node B               see Rossier Equation (2.22)
    end
    
end

function [RF,RM,L,T,Ma,Mb] = DrResiduals(RF,Links,Type,X,Young,A,Iy,L0,A_nodes,A_a0,A_b0,P)

% Definition of local varibales
nElements = size(Links,1);
nDim = size(Links,2)/2;                  % problem dimension
L = L0*0;                            % Initialize L
AElement = zeros(nElements,1);
A_a = A_a0*0;                        % Initialize Angles
A_b = A_b0*0;
RM = zeros(nElements+1,1);

if nDim == 3
    for i = 1:nElements		% for all links
        L(i) = sqrt((X(Links(i,4))-X(Links(i,1)))^2 + (X(Links(i,5))-X(Links(i,2)))^2 + (X(Links(i,6))-X(Links(i,3)))^2);
    end
else
    for i = 1:nElements		% for all links
        delta = [X(Links(i,3))-X(Links(i,1)) X(Links(i,4))-X(Links(i,2))];
        AElement(i) = atan(delta(2)/delta(1));
        A_a(i) = A_nodes(i) - AElement(i);
        A_b(i) = A_nodes(i+1) - AElement(i);
        L(i) = sqrt(delta(1)^2 + delta(2)^2);
    end    
end

% Compute internal forces
T = ((Young.*A)./L0).*(L-L0) + P;

Ma = (2*Young.*Iy./L).*(2*A_a + A_b);
Mb = (2*Young.*Iy./L).*(2*A_b + A_a);

         
Va = (Ma-Mb)./L;

% Remove compression in cables
T = T.*~(Type.*(T < 0));


for i = 2:nElements
    
    RM(i) = Ma(i)+Mb(i-1);
    
end

for i = 1:nElements
    % residual force computation
    % for all direction
    
    for j = 1:nDim
        rF = (T(i)/L(i))*(X(Links(i,j+nDim))-X(Links(i,j)));	% residual force = (internal force / link length) * Delta coordinate see Barnes Equation (4) / Rossier Equation (2.17)
        
        if j == 1
        rV = Va(i)*sin(AElement(i));
        else
        rV = Va(i)*cos(AElement(i));   
        end
        
        RF(Links(i,j)) = RF(Links(i,j))+rF+rV;			% residual force = external force (see above) + residual force                                                                      see Barnes Equation (4) / Rossier Equation (2.19)
        RF(Links(i,j+nDim)) = RF(Links(i,j+nDim))-rF-rV;			% residual force = external force (see above) - residual force
        % in one node is r added and in the other r is substracted
    end

end