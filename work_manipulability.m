function s = work_manipulability(obj,samplesNo)

% WORKMANIPULABILITY calculates the manipulability of the robot for
% sampleNo randomly generated poses 
%
% SYNTAX s = work_manipulability(obj,samplesNo)
%
% INPUT obj: robot object
% samplesNo: number of sample poses
%   
% OUTPUT  s: structure containing position, orientation and 
% manipulability measure for each entry
%
% EXAMPLES: m = work_manipulability(robot,3000)
%
% TO DO:
% - preallocate variables
% - implement extended measure.

dens = 100;
l{1} = linspace(0, 0, dens);
l{2} = linspace(0, 0, dens);
l{3} = linspace(0, 0, dens);
l{4} = linspace(-170*3.14/180, 170*3.14/180, dens);
l{5} = linspace(-120*3.14/180, 120*3.14/180, dens);
l{6} = linspace(-170*3.14/180, 170*3.14/180, dens);
l{7} = linspace(-120*3.14/180, 120*3.14/180, dens);
l{8} = linspace(-170*3.14/180, 170*3.14/180, dens);
l{9} = linspace(-120*3.14/180, 120*3.14/180, dens);
l{10} = linspace(-170*3.14/180, 170*3.14/180, dens);

i=1;
curAngles = zeros(1,size(obj.links,2)-3);

while i<samplesNo
    for j = 1:size(obj.links,2)-1
        curAngles(j) = l{j}(randi(dens));
    end
    T = fkNum(curAngles);
    J = jacobianNum(curAngles);
    J = J(:,4:end);
    s(i) = struct('pos', T(1:3,4), 'ori', T(1:3,1:3), ...
        'manip', manipulability(obj,J,curAngles,1));
    i=i+1;
end


for i = 1:size(s,2)
   x(:,i) = s(i).pos;
end
m = [s(:).manip];
m_max = max(m);
m_min = min(m);

% scatter3(x(1,x(2,:)>0.16),x(2,x(2,:)>0.16),x(3,x(2,:)>0.16),100,m(x(2,:)>0.16)./(m_max),'.');
scatter3(x(1,:),x(2,:),x(3,:),100,m(:)./(m_max),'.'); %plot all

end