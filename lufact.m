function lufact(A,b)
% Solve the system Ax=b using the LU decomposition.
n=length(b);
y=zeros(n,1);
x=zeros(n,1);
fprintf('\n');

U=eye(n);

L(:,1) = A(:,1)/U(1,1);
U(1,:) = A(1,:)/L(1,1);

for i=2:n-1
    S = L(i,1:i-1) * U(1:i-1,i);
    L(i,i)=(A(i,i)-S)/U(i,i);

    for j=i+1:n
        S = L(j,1:i-1) * U(1:i-1,i);
        L(j,i)=(A(j,i)-S)/U(i,i);

        S = L(i,1:i-1) * U(1:i-1,j);
        U(i,j)=(A(i,j)-S)/L(i,i);
    end
end

S = L(n,1:n-1) * U(1:n-1,n);
L(n,n)=(A(n,n)-S)/U(n,n);
% Perform the forward substitution.
y(1)=b(1)/L(1,1);
for i=2:n
    S = b(i) - L(i,1:i-1) * y(1:i-1);
    y(i) = S / L(i,i);
end
% Perform the back substitution.
x(n)=y(n)/U(n,n);
for i=n-1:-1:1
    S = y(i) - U(i,i+1:n) * x(i+1:n);
    x(i)=S/U(i,i);
end
% Print the results
L
disp('         The forward substitution gives')
y
U
disp('         The vector solution is =')
x
