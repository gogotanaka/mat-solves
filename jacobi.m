function jacobi(A,b,x0,tol,itmax)
% Solve the system Ax=b using Jacobi iteration method.
n=length(b);
x=zeros(n,1);
fprintf('\n');
disp('      The augmented matrix is =')
Augm=[A b]
Y=zeros(n,1);
Y=x0;
for k=1:itmax+1
    U = triu(A,1);
    D = diag(diag(A));
    L = tril(A,-1);

    if ismember(0,diag(D))
        break
    end

    x = D^(-1) * (b - (L + U) * x0);

    err=abs(norm(x-x0));
    rerr=err/(norm(x)+eps);
    x0=x;
    Y=[Y x];
    if(rerr<tol)
        break
    end
end



% Print the results
if ismember(0,diag(D))
    disp('        division by zero')
elseif (k==itmax+1)
    disp('        No convergence')
else
    fprintf('\n');
    disp(' The solution vectors are:')
    fprintf('\n');
    disp('iter #       0         1         2         3        4        ...')
    fprintf('\n');
    for i=1:n
        fprintf('x%1.0f =    ',i)
        fprintf('%10.6f ',Y(i,[1:k+1]))
        fprintf('\n');
    end
    fprintf('\n');
    disp(['The method converges after ',num2str(k),' iterations to']);
    x
end
