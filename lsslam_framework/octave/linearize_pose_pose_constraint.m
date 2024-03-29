% Compute the error of a pose-pose constraint
% x1 3x1 vector (x,y,theta) of the first robot pose
% x2 3x1 vector (x,y,theta) of the second robot pose
% z 3x1 vector (x,y,theta) of the measurement
%
% You may use the functions v2t() and t2v() to compute
% a Homogeneous matrix out of a (x, y, theta) vector
% for computing the error.
%
% Output
% e 3x1 error of the constraint
% A 3x3 Jacobian wrt x1
% B 3x3 Jacobian wrt x2
function [e, A, B] = linearize_pose_pose_constraint(x1, x2, z)

  % TODO compute the error and the Jacobians of the error
  t1 = v2t(x1);
  t2 = v2t(x2);
  t12 = v2t(z);
  e = t2v(invt(t12) * (invt(t1) * t2));
  B = diag([0,0,1]);
  B(1:2,1:2) = t12(1:2,1:2)' * t1(1:2,1:2)';
  A = -B;
  A(1:2,3) = -B(1:2,1:2) * [0 -1;1 0] * (x2(1:2) - x1(1:2));
end;
