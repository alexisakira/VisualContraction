function [Cmat,Vmat] = get_cv(sg,Vmat0)
% update optimal consumption and value function from Bellman equation

aGrid = sg.aGrid;
N = length(aGrid); % number of grid points
nz = size(sg.P,1); % number of states

Cmat = zeros(nz,N);
Vmat = zeros(nz,N);

for z = 1:nz
    for n = 1:N
        a = aGrid(n);
        if a == 0
            c = 0;
            Cmat(z,n) = c;
            Vmat(z,n) = get_cv_obj(c,a,z,sg,Vmat0);
        else
            func = @(c)(-get_cv_obj(c,a,z,sg,Vmat0));
            [c,fval] = fminbnd(func,0,a); % compute optimal consumption
            Cmat(z,n) = c;
            Vmat(z,n) = -fval;
        end
    end
end

end

