function w = manipulability_ext(obj,Jt,q,use_svd)

% MANIPULABILITY calculates the manipulability of the robot for a given pose
%
% SYNTAX w = manipulability(obj,Jt,q,use_svd)
%
% INPUT obj: robot object
%        Jt: robot numeric jacobian
%         q: vector of the joint angles in radiants
%   use_svd: use svd decomposition, otherwise use sqrt(det(JJ'))
%
% OUTPUT  w: manipulability measure
%
% EXAMPLES: m = manipulability(robot,[1,2,3,...,0],[0,0,0,...,0],1)
%
% TO DO:
% -jacobian sub-matrix selection hardcoded
%   change to w = sqrt(det(j*j')) to use full jacobian.

selfCollision = 0;
jointLimits = 0;
collision = 1;
PLOT = 0;
curAngles = q;
J = Jt(:,4:end);
T = fkNum([0,0,0,curAngles]);  % forward kinematics
L = ones(size(J,1),size(J,2)); % joint limits penalization matrix
O = ones(size(J,1),size(J,2)); % self collision penalization matrix

dir = [0,1,0,0,0,0];
gamma = dir/norm(dir,2);

l{1} = [-170*3.14/180, 170*3.14/180];
l{2} = [-120*3.14/180, 120*3.14/180];
l{3} = l{1};
l{4} = l{2};
l{5} = l{1};
l{6} = l{2};
l{7} = l{1};

if (jointLimits)
    p_minus = ones(1,7);
    p_plus = ones(1,7);
    for j = 1:7
        if (abs(curAngles(j) - l{j}(1)) > abs(curAngles(j) - l{j}(end)))
            p_minus(j) = 1;
            p_plus(j) = 1/sqrt(1+abs(grad_jnt_lim(l{j}(end),l{j}(1),curAngles(j))));
        else
            p_minus(j) = 1/sqrt(1+abs(grad_jnt_lim(l{j}(end),l{j}(1),curAngles(j))));
            p_plus(j) = 1;
        end
    end
    for i = 1:size(J,1)
        for j = 1:size(J,2)
            if (sign(J(i,j))*sign(gamma(i)) < 0)
                L(i,j) = p_minus(j);
            elseif (sign(J(i,j))*sign(gamma(i)) > 0)
                L(i,j) = p_plus(j);
            else
                L(i,j) = 1;
            end
        end
    end
end

J = Jt;
shelf = loadshelf([0 -0.0]);

if (collision)
    convmeshes = updateRobot(obj, [0,0,0,curAngles], obj.convmeshes);
    [coll,dist,Pa,Pb] = gjk(shelf,convmeshes(11),12);
    dist
    T1 = fk11([0,0,0,curAngles]);
    T2(1:3,1:3) = T1(1:3,1:3);
    T2(1:3,4) = Pa;
    T2(4,1:4) = [0,0,0,1];
    Ja = zeros(6,10);
    Jb = partialJacobian(obj,[0,0,0,curAngles],T2,11);
    pa = [Pa';0;0;0];
    pb = [Pb';0;0;0];   
    % grad_distance(alpha, beta, rho, pa, pb, Ja, Jb)
    grad = grad_distance(50, 1/2, 30, pa, pb, Ja, Jb);
    if PLOT
        convmeshes = updateRobot(obj, [0,0,0,curAngles(1,:)], obj.convmeshes);
        meshes = updateRobot(obj, [0,0,0,curAngles(1,:)], obj.meshes);
        obj.plotRobot(meshes,convmeshes,0.5);
        line([pa(1) pb(1)],[pa(2) pb(2)],[pa(3) pb(3)],'LineWidth',2)
        %pause
        drawnow;
    end
        
    v = pa - pb;
    v = [v;0;0;0];

    J = Jt(:,4:end);

    % construct O
    for i = 1:size(J,1)
        for j = 1:size(J,2)
            if (sign(gamma(i)) < 0)
                if (v(i) > 0 && i > 3)
                    o_minus = 1;
                else
                    o_minus = 1/sqrt(1+abs(grad(j+3)));
                end
                O(i,j) = o_minus;
            elseif (sign(gamma(i)) > 0)
                if (v(i) > 0 && i <= 3)
                    o_plus = 1/sqrt(1+abs(grad(j+3)));
                else
                    o_plus = 1;
                end
                O(i,j) = o_plus;
            else
                O(i,j) = 1;
            end
        end
    end
end

if (selfCollision)
    collMatrTest = triu(ones(size(obj.links,2)),2);
    collMatrTest(1:2,:) = 0;
    collMatrTest(:,1:2) = 0;
    collMatrTest(:,4) = 0;
    collMatrTest(9,11) = 0;
    %     collMatrTest = zeros(size(obj.links,2));
    %     collMatrTest(3,11) = 1; %only base=3, end effector=11
    
    convmeshes = updateRobot(obj, [0,0,0,curAngles], obj.convmeshes); % update
    % meshes = updateRobot(obj, curAngles, obj.meshes); % update
    %     mesh3.vertices = convmeshes(3).v;
    %     mesh3.faces = convmeshes(3).f;
    %     mesh11.vertices = convmeshes(11).v;
    %     mesh11.faces = convmeshes(11).f;
    
    flag = 0;
    for i = 1:size(obj.links,2) % GJK algorithm
        for j = i:size(obj.links,2)
            if collMatrTest(i,j) == 1
                [coll,dist(i,j),Pa{i,j},Pb{i,j}] = gjk(convmeshes(i),convmeshes(j),12);
                if coll
                    flag = 1;
                    break;
                end
            end
        end
    end % GJK algorithm
    
    %     [coll,dist(3,11),Pa{3,11},Pb{3,11}] = gjk(convmeshes(3),convmeshes(11),6);
    
    if (~flag) % no links in collision state
        %         row = 3;
        %         col = 11;
        %construct O
        [row,col] = find(dist == min(dist(dist>0)));
        pa = Pa{row,col}'; % point on the shape i (mobile base)
        T1 = feval(strcat('fk',int2str(row)),([0,0,0,curAngles(1:row-4)])); %hardcoded
        T2(1:3,1:3) = T1(1:3,1:3);
        T2(1:3,4) = pa;
        T2(4,1:4) = [0,0,0,1];
        Ja = partialJacobian(obj,[0,0,0,curAngles],T2,row);
        %         Ja = Jt(:,1:row-1);
        %         Ja = obj.kinematicTransformation(Ja,T1,T2);
        
        pb = Pb{row,col}'; % point on the shape j (end effector)
        %             pb = Pb'; % point on the shape j (end effector)
        T1 = feval(strcat('fk',int2str(col)),([0,0,0,curAngles(1:col-4)])); %hardcoded
        T2(1:3,1:3) = T1(1:3,1:3);
        T2(1:3,4) = pb;
        T2(4,1:4) = [0,0,0,1];
        Jb = partialJacobian(obj,[0,0,0,curAngles],T2,col);
        %         Jb = Jt(:,1:col-1);
        %         Jb = obj.kinematicTransformation(Jb,T1,T2);
        
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
        % grad_distance(alpha, beta, rho, pa, pb, Ja, Jb)
        grad = grad_distance(50, 1/2, 30, pa, pb, Ja, Jb);
        
        if PLOT
            convmeshes = updateRobot(obj, [0,0,0,curAngles(1,:)], obj.convmeshes);
            meshes = updateRobot(obj, [0,0,0,curAngles(1,:)], obj.meshes);
            obj.plotRobot(meshes,convmeshes,0.5);
            line([pa(1) pb(1)],[pa(2) pb(2)],[pa(3) pb(3)],'LineWidth',2)
            %pause
            drawnow;
        end
        
        v = pa - pb;
        v = [v;0;0;0];
        
        J = Jt(:,4:end);
        
        % construct O
        for i = 1:size(J,1)
            for j = 1:size(J,2)
                if (sign(gamma(i)) < 0)
                    if (v(i) > 0 && i > 3)
                        o_minus = 1;
                    else
                        o_minus = 1/sqrt(1+abs(grad(j+3)));
                    end
                    O(i,j) = o_minus;
                elseif (sign(gamma(i)) > 0)
                    if (v(i) > 0 && i <= 3)
                        o_plus = 1/sqrt(1+abs(grad(j+3)));
                    else
                        o_plus = 1;
                    end
                    O(i,j) = o_plus;
                else
                    O(i,j) = 1;
                end
            end
        end
    end
end

J = Jt(:,4:end);

J = L.*O.*J;
w = manipulability(obj,J,curAngles,use_svd);

end