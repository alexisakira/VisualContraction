function os = solve_os(os)
%% solve optimal savings problem by policy function iteration

nz = size(os.P,1); % number of states
N = length(os.aGrid); % number of grid points
MaxIter = os.MaxIter; % maximum number of iterations
tol = os.tol; % error tolerance

CMat = zeros(nz,N,MaxIter); % store consumption function for plotting later
Cmat0 = os.Cmat0; % initialize consumption function
CMat(:,:,1) = Cmat0; 

% policy function iteration
for i = 1:MaxIter
    Cmat = get_c(os,Cmat0);
    CMat(:,:,i+1) = Cmat; 
    if max(max(abs(Cmat./Cmat0 - 1))) < tol
        str = 'Converged after %3.0f iterations\n';
        fprintf(str,i)
        os.Cmat = Cmat; % consumption function
        os.CMat = CMat(:,:,1:i+1); % remove irrelevant part of VMat
        break
    else % update consumption function
        Cmat0 = Cmat;
    end
end

os.imax = i; % number of iterations required for convergence

end

