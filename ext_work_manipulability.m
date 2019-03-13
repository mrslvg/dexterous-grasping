function s = ext_work_manipulability(obj,samplesNo,jointLimits,selfCollision)

% EXTWORKMANIP calculates the extended manipulability of the robot for
% sampleNo randomly generated poses
%
% SYNTAX s = ext_work_manip(obj,samplesNo)
%
% INPUT       obj: robot object
%       samplesNo: number of sample poses
%     jointLimits: perform joint limits penalization
%   selfCollision: perform self-collision penalization
%
% OUTPUT        s: structure containing position, orientation and
% manipulability measure for each entry
%
% EXAMPLES: m = ext_work_manipulability(robot,3000,true,true)
%
% TO DO:
% - preallocate variables

A = importdata('RobotData.txt', ' ');

PLOT = 0;

if (selfCollision)
    collMatrTest = triu(ones(size(obj.links,2)),2);
    collMatrTest(1:2,:) = 0;
    collMatrTest(:,1:2) = 0;
    collMatrTest(4,:) = 0;
    collMatrTest(:,4) = 0;
    collMatrTest(9,11) = 0;
    %     collMatrTest = zeros(size(obj.links,2));
    %     collMatrTest(3,11) = 1; %only base=3, end effector=11
end

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

% initialize variables
iterations=1;
curAngles = zeros(1,size(obj.links,2)-1);
dir = [-1,1,-1,-1,-1,-1];
gamma = dir/norm(dir,2);
if jointLimits
    p_minus = ones(1,size(obj.links,2)-1);
    p_plus = ones(1,size(obj.links,2)-1);
end

% while iterations<samplesNo
while iterations<size(A.data,1)-100
    fprintf('iteration: %s\n',int2str(iterations));
    
    for j = 4:size(obj.links,2)-1
%         curAngles(j) = l{j}(randi(dens)); % generate random angles
        curAngles(j) = A.data(iterations,j-2); % generate random angles

        if jointLimits
            if (abs(curAngles(j) - l{j}(1)) > abs(curAngles(j) - l{j}(end)))
                p_minus(j) = 1;
                p_plus(j) = 1/sqrt(1+abs(grad_jnt_lim(l{j}(end),l{j}(1),curAngles(j))));
            else
                p_minus(j) = 1/sqrt(1+abs(grad_jnt_lim(l{j}(end),l{j}(1),curAngles(j))));
                p_plus(j) = 1;
            end
        end
    end
    
    T = fkNum(curAngles); % forward kinematics
    
    convmeshes = updateRobot(obj, curAngles, obj.convmeshes); % update
%     meshes = updateRobot(obj, curAngles, obj.meshes); % update
    mesh3.vertices = convmeshes(3).v;
    mesh3.faces = convmeshes(3).f;
    mesh11.vertices = convmeshes(11).v;
    mesh11.faces = convmeshes(11).f;
    
    in = inpolyhedron(mesh3, mesh11.vertices);
    
    if (~in)
        
        J = jacobianNum(curAngles); % Jacobian
        L = ones(size(J,1),size(J,2)-3); % joint limits penalization matrix
        O = ones(size(J,1),size(J,2)-3); % self collision penalization matrix
        
        if selfCollision % self collision
            flag = 0;
%             for i = 1:size(obj.links,2) % GJK algorithm
%                 for j = i:size(obj.links,2)
%                     if collMatrTest(i,j) == 1
%                         [coll,dist(i,j),Pa{i,j},Pb{i,j}] = gjk(convmeshes(i),convmeshes(j),12);
%                         if coll
%                             flag = 1;
%                             break;
%                         end
%                     end
%                 end
%             end % GJK algorithm
            
            [coll,dist(3,11),Pa{3,11},Pb{3,11}] = gjk(convmeshes(3),convmeshes(11),6);

            
            if (~flag) % no links in collision state
                row = 3;
                col = 11;
                %construct O
                %[row,col] = find(dist == min(dist(dist>0)));
                pa = Pa{row,col}'; % point on the shape i (mobile base)
                Ja = J(:,1:row);
%                 T1 = fk3(curAngles); %hardcoded
%                 T2(1:3,1:3) = eye(3);
%                 T2(1:3,4) = pa;
%                 T2(4,1:4) = [0,0,0,1];
%                 Ja = obj.kinematicTransformation(Ja,T1,T2);
                
                pb = Pb{row,col}'; % point on the shape j (end effector)
                %             pb = Pb'; % point on the shape j (end effector)
                Jb = J(:,1:col-1);
                T1 = T; %hardcoded
                T2(1:3,1:3) = eye(3);
                T2(1:3,4) = pb;
                T2(4,1:4) = [0,0,0,1];
                Jb = obj.kinematicTransformation(Jb,T1,T2);
                
                %             Ja = double(vpa(subs(Ja, ...
                %                 {'q1','q2','q3','q4','q5','q6','q7','q8','q9','q10'}, ...
                %                 {curAngles(1),curAngles(2),curAngles(3),curAngles(4), ...
                %                 curAngles(5),curAngles(6),curAngles(7),curAngles(8), ...
                %                 curAngles(9),curAngles(10)})));
                
                %             Jb = double(vpa(subs(Jb, ...
                %                 {'q1','q2','q3','q4','q5','q6','q7','q8','q9','q10'}, ...
                %                 {curAngles(1),curAngles(2),curAngles(3),curAngles(4), ...
                %                 curAngles(5),curAngles(6),curAngles(7),curAngles(8), ...
                %                 curAngles(9),curAngles(10)})));
                
                Ja = [Ja,zeros(6,size(obj.links,2)-1-size(Ja,2))];
                Jb = [Jb,zeros(6,size(obj.links,2)-1-size(Jb,2))];
                
                pa = [pa;0;0;0];
                pb = [pb;0;0;0];
                grad = grad_distance(1, 1, 0.2, pa, pb, Ja, Jb);
                
                if PLOT
                    obj.plotRobot(meshes,convmeshes);
                    line([pa(1) pb(1)],[pa(2) pb(2)],[pa(3) pb(3)],'LineWidth',2)
                    %pause
                    drawnow;
                end
                
                v = pa - pb;
                v = [v;0;0;0];
                
                J = J(:,4:end);
                
                % construct O
                for i = 1:size(J,1)
                    for j = 1:size(J,2)
                        if (sign(gamma(i)) < 0)
                            if (v(i) > 0 || i > 3)
                                o_minus = 1;
                            else
                                o_minus = 1/sqrt(1+abs(grad(j+3)));
                            end
                            O(i,j) = o_minus;
                        else
                            if (v(i) > 0 && i <= 3)
                                o_plus = 1/sqrt(1+abs(grad(j+3)));
                            else
                                o_plus = 1;
                            end
                            O(i,j) = o_plus;
                        end
                    end
                end
                
                J = L.*O.*J;
                
                s(iterations) = struct('pos', T(1:3,4), 'ori', T(1:3,1:3), ...
                    'manip', manipulability(obj,J,curAngles,1));
                iterations=iterations+1;
                
            end % no links in collision state
        end % self collision
        
        if jointLimits % joint limits
            % construct L
            for i = 1:size(J,1)
                for j = 1:size(J,2)
                    if (sign(J(i,j))*sign(gamma(i)) < 0)
                        L(i,j) = p_minus(j);
                    else
                        L(i,j) = p_plus(j);
                    end
                end
            end
        end % joint limits
        
    end % end effector not in collision with base
    
end

for i = 1:size(s,2)
    x(:,i) = s(i).pos;
end
m = [s(:).manip];
m_max = max(m);
m_min = min(m);
% scatter3(x(1,x(2,:)>0.16),x(2,x(2,:)>0.16),x(3,x(2,:)>0.16),100,m(x(2,:)>0.16)./(m_max),'.');
%scatter3(x(1,:),x(2,:),x(3,:),100,m(:)./(m_max),'.'); %plot all

end