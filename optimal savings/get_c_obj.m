function obj = get_c_obj(c,a,z,os,Cmat0)
% objective function for computing optimal consumption
% c: current consumption
% a: current asset
% z: current state
% os: optimal savings structure
% cmat0: matrix that stores comsumption functions

nz = size(os.P,1); % number of states
obj = os.mu(c,z); % marginal utility

% compute difference of Euler equation
for z1 = 1:nz
    R1 = os.R(z,z1); % gross return
    a1 = R1*(a-c) + os.Y(z,z1); % next period's asset
    c1 = interp1(os.aGrid,Cmat0(z1,:),a1,'linear','extrap'); % next period's consumption
    obj = obj - os.P(z,z1)*os.beta(z,z1)*R1*os.mu(c1,z1);
end

end

