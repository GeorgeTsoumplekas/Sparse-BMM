clc;
clear;

n=5e4;
d=2;
A=sprand(n,n,d/n)>0;
B=sprand(n,n,d/n)>0;
tic;
C=(A*B)>0;
toc

% A=[0,0,0,0;
%     0,0,1,0;
%     1,0,0,1;
%     1,0,0,0];
% 
% B=[0,1,0,1;
%     0,0,0,0;
%     0,0,1,0;
%     0,1,0,0;];

tic;

%Extract A in COO format
nnz_A = nnz(A);
[n,~] = size(A);
fileID = fopen('A_coo.txt','w');
fprintf(fileID,"%d,%d\n",n,nnz_A);
[row, col, ~] = find(A');
for i=1:nnz_A
   fprintf(fileID,"%d,%d\n",col(i),row(i)); 
end

%Extract B in COO format
nnz_B = nnz(B);
[n,~] = size(B);
fileID = fopen('B_coo.txt','w');
fprintf(fileID,"%d,%d\n",n,nnz_B);
[row, col, ~] = find(B);
for i=1:nnz_B
   fprintf(fileID,"%d,%d\n",row(i),col(i)); 
end

%Extract C in COO format
nnz_C = nnz(C);
[n,~] = size(C);
fileID = fopen('C_coo.txt','w');
fprintf(fileID,"%d,%d\n",n,nnz_C);
[row, col, ~] = find(C');
for i=1:nnz_C
   fprintf(fileID,"%d,%d\n",col(i),row(i)); 
end

toc