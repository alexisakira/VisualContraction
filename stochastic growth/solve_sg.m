function sg = solve_sg(sg)
%% solve stochastic growth model by value function iteration

nz = size(sg.P,1); % number of states
N = length(sg.aGrid); % number of grid points
MaxIter = sg.MaxIter; % maximum number of iterations
tol = sg.tol; % error tolerance
nu = sg.nu; % parameter to skip maximization step

VMat = zeros(nz,N,MaxIter); % store value function for plotting later
Vmat0 = sg.Vmat0; % initialize value function
VMat(:,:,1) = Vmat0; 

% value function iteration
for i = 1:MaxIter
    if rem(i-1,nu) == 0 % carry out maximization
        [Cmat,Vmat] = get_cv(sg,Vmat0);
    else % no maximization
        for z = 1:nz
            for n = 1:N
                a = sg.aGrid(n);
                c = Cmat0(z,n);
                Vmat(z,n) = get_cv_obj(c,a,z,sg,Vmat0);
            end
        end
    end
    VMat(:,:,i+1) = Vmat; 
    if max(max(abs(Vmat - Vmat0))) < tol
        str = 'Converged after %3.0f iterations\n';
        fprintf(str,i)
        sg.Cmat = Cmat; % consumption function
        sg.Vmat = Vmat; % value function
        sg.VMat = VMat(:,:,1:i+1); % remove irrelevant part of VMat
        break
    else % update value function
        Cmat0 = Cmat;
        Vmat0 = Vmat;
    end
end

sg.imax = i; % number of iterations required for convergence

end

