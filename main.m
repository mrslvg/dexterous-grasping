close all
clear
clc

joints = importdata('joints.txt');
joints(:,4) = -joints(:,4);
platform = importdata('platform.txt');

i = 1;
while i<=size(joints,1)
    T = fkNum([0,0,0,joints(i,:)]);
    J = jacobianNum([0,0,0,joints(i,:)]);
    J = J(:,4:end);
    s(i) = struct('pos', platform(i,:), ...
        'manip', manipulability(r,J,joints(i,:),1));
    s_ext(i) = struct('pos', platform(i,:), ...
        'manip', manipulability_ext(r,J,joints(i,:),1));
    i=i+1;
end