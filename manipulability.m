function w = manipulability(obj,J,q,use_svd)

% MANIPULABILITY calculates the manipulability of the robot for a given pose
%
% SYNTAX w = manipulability(obj,q,use_svd)
%
% INPUT obj: robot object
%         J: robot numeric jacobian
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

jacob = J;%jacobianNum(q);

if use_svd
    
    manip_measure = 1;
    
    try
        singular_values = svd(jacob*jacob');
        for i = 1 : size(singular_values,1)
            manip_measure  = manip_measure * singular_values(i);
        end
        w = min(singular_values)/max(singular_values);
    catch
        w = 0;
    end
    
else
    
    w = sqrt(det(jacob*jacob'));
    
end
end