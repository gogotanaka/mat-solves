A = [1 1 1 1; 2 3 1 5; -1 1 -5 3; 3 1 7 -2];
b = [10 31 -2 18];

lufact(A,b)


A = [7 -2 1 0; 1 -9 3 -1; 2 0 10 1; 1 -1 1 6];
b = [17 13 15 10]';
x0 = [0 0 0 0];
jacobi(A, b, x0, 10^(-3), 30)