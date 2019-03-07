function X=rownorm(X)
% normalize rows to unit length.
s=abs(sqrt(sum(X.^2,2)));

if any(s==0),
  warning('Contains zero vectors: can''t normalize them!');
end

X=X./repmat(s,1,size(X,2));