function X = normalize_columns(X)
    X=abs(X);
    col_sum = sum(X, 1);
    col_sum(col_sum == 0) = 1; 
    X = X ./ col_sum; 
end