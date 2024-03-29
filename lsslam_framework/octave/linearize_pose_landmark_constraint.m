% Compute the error of a pose-landmark constraint
% x 3x1 vector (x,y,theta) of the robot pose
% l 2x1 vector (x,y) of the landmark
% z 2x1 vector (x,y) of the measurement, the position of the landmark in
%   the coordinate frame of the robot given by the vector x
%
% Output
% e 2x1 error of the constraint
% A 2x3 Jacobian wrt x
% B 2x2 Jacobian wrt l
function [e, A, B] = linearize_pose_landmark_constraint(x, l, z)

  % TODO compute the error and the Jacobians of the error
  t = v2t(x);
  Rx = t(1:2,1:2);
  e = Rx' * (l - x(1:2)) - z; 
  B = Rx';
  A = zeros(2,3);
  A(1:2,1:2) = -Rx';
  A(1:2,3) = -Rx' * [0 -1;1 0] * (l - x(1:2));
end;
