% initialize model
install; 
Assembly = BaseLayout(modelDataFolder);

% assign load
Assembly.ExternalLoad(2:2:22) = -10;

% run DR
[Xall,Tall,Lall,N_nodes,Ma,Mb] = DR(Assembly);

Xall = Xall(:,end);

%% plot results
figure;
for i = 1:Assembly.nElements
    plot([Assembly.Coordinates(i,1) Assembly.Coordinates(i+1,1)],[Assembly.Coordinates(i,2) Assembly.Coordinates(i+1,2)],'-k')
    hold on
end

CoordinatesDeformed = reshape(Xall,Assembly.nDim,Assembly.nNodes)';

for i = 1:Assembly.nElements
    plot([CoordinatesDeformed(i,1) CoordinatesDeformed(i+1,1)],[CoordinatesDeformed(i,2) CoordinatesDeformed(i+1,2)],'-r')
    hold on
end

axis equal