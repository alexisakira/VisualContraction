function Cmat = get_c(os,Cmat0)
% update optimal consumption from Euler equation

aGrid = os.aGrid;
N = length(aGrid); % number of grid points
nz = size(os.P,1); % number of states

Cmat = zeros(nz,N);

for z = 1:nz
    for n = 1:N
        a = aGrid(n);
        if a == 0
            Cmat(z,n) = 0;
        else % a > 0
            func = @(c)(get_c_obj(c,a,z,os,Cmat0)); % difference of Euler equation
            if func(0) <= 0
                Cmat(z,n) = 0; % consuming nothing is optimal
            elseif func(a) >= 0
                Cmat(z,n) = a; % consuming everything is optimal
            else
                x0 = a*[1e-3 1];
                c = fzero(func,x0); % compute optimal consumption
                Cmat(z,n) = c;
            end
        end
    end
end

end

