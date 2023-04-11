function T= trophiclevel(A)
%T= trophiclevel(A)
%	Computes trophic level of each species in adjacency matrix A

websize= size(A,1);

% Remove self-loops
for j= 1:websize A(j,j)= 0; end

% Number of prey each species has (sum of columns)
nprey= sum(A,1);

% Unweighted Q matrix.
Q= zeros(websize);
[row,col]= find(A~=0);
for j= 1:length(row)
	Q(row(j),col(j))= A(row(j),col(j))/nprey(col(j));
end

% Calculate trophic levels as T2=(I-Q)^-1 * 1 (Levine 1980 JTB Vol 83) (If the matrix
% is not singular)
% if rcond(A) < eps
% 	T= NaN(websize,1);
% else
	T=(inv(eye(websize)-Q'))*ones(websize,1);
% end
