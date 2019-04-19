function [L, P] = lqr_infinite_horizon_solution(Q, R)

%% find the infinite horizon L and P through running LQR back-ups
%%   until norm(L_new - L_current, 2) <= 1e-4  
dt = 0.1;
mc = 10; mp = 2.; l = 1.; g= 9.81;

% TODO write A,B matrices
a1 = mp * g / mc; 
a2 = (mc + mp) * g / l / mc;

A = eye(4) + dt * [0 0 1 0; 0 0 0 1; 0 a1 0 0; 0 a2 0 0];
B = dt * [0; 0; 1 / mc; 1/ l / mc];

% TODO implement Riccati recursion
Q = eye(4);
R = eye(1);
[P1,L1,G1] = dare(A,B,Q,R)
K = lqr (A,B,Q,R)
P = Q;
tranA = transpose(A);
tranB = transpose(B);
G = inv(tranB*P*B + R)*tranB*P*A; 
L_current = eig(A-B*G);
L_new = L_current;

while  1
	L_current = L_new;
	P = Q + tranA*P*A - tranA*P*B * G;
	G = inv(tranB*P*B + R)*tranB*P*A;
	L_new = eig(A-B*G);
	if (norm(L_new - L_current, 2) <= 1e-4) 
		break;
	end
end

return