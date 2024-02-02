function Out = get_cv_obj(c,a,z,sg,Vmat0)
% objective function for computing optimal consumption
% c: current consumption
% a: current resource
% z: current state
% sg: stochastic growth structure
% Vmat0: matrix that stores value functions

nz = size(sg.P,1); % number of states
V = sg.u(c,z); % flow utility

for z1 = 1:nz
    a1 = sg.f(a-c,z1); % next period's resource
    V = V + sg.beta*sg.P(z,z1)*interp1(sg.aGrid,Vmat0(z1,:),a1,'spline');
end

Out = V;

end

