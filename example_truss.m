% initialize model
install; 
Assembly = BaseLayout(modelDataFolder);

% assign load

load =-2;

for i=1:Assembly.nNodes
    if Assembly.Coordinates(i,2)==100
        Assembly.ExternalLoad(2*i)=load;
    end
end

% run DR
[Xall,Tall,Lall] = DR(Assembly);

%% plot results
figure;
for i = 1:Assembly.nElements
    n1=Assembly.Elements2Nodes(i,1);
    n2=Assembly.Elements2Nodes(i,2);
    plot([Assembly.Coordinates(n1,1) Assembly.Coordinates(n2,1)],[Assembly.Coordinates(n1,2) Assembly.Coordinates(n2,2)],'-k')
    hold on
end

CoordinatesDeformed = reshape(Xall(:,end),Assembly.nDim,Assembly.nNodes)';

for i = 1:Assembly.nElements
    n1=Assembly.Elements2Nodes(i,1);
    n2=Assembly.Elements2Nodes(i,2);
    plot([CoordinatesDeformed(n1,1) CoordinatesDeformed(n2,1)],[CoordinatesDeformed(n1,2) CoordinatesDeformed(n2,2)],'-r')
    hold on
end

% % Animation de l'onde
%  filename = 'dep.gif';
%  for n = 1:127
%       CoordinatesDeformed = reshape(Xall(:,i),Assembly.nDim,Assembly.nNodes)';
%       for i = 1:Assembly.nElements
%       n1=Assembly.Elements2Nodes(i,1);
%       n2=Assembly.Elements2Nodes(i,2);
%       plot([CoordinatesDeformed(n1,1) CoordinatesDeformed(n2,1)],[CoordinatesDeformed(n1,2) CoordinatesDeformed(n2,2)],'-r')
%       hold on
%       end
%       %scatter(1:nb_noeuds,U(:,n),'filled')
%       %axis ([0 nb_noeuds min(min(U)) max(max(U))])
%       drawnow
%       frame = getframe(1);
%       im = frame2im(frame);
%       [imind,cm] = rgb2ind(im,256);
%       if n == 1;
%           imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1/240);
%       else
%           imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1/240);
%       end
%  end

axis equal