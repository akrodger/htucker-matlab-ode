function [mse_val] = ht_mse(X,Y)
%
%    Mean Square Error of HT tensors by Bram Rodgers
%    Original Draft 28 Oct, 2021
%
%
%    Description of This Function:
%        Computes the mean square error metric between two HTucker
%        factorizations by using a the formula
%        
%        norm(X-Y)/sqrt(numel(X))
%
%	Argument List:
%        X:  A HTucker factored tensor
%        Y:  A HTucker factored tensor
%
%    Return List:
%        mse_val: The mean square error between X and Y
%
    sz = size(X);
    num_elements = prod(sz);
    D = minus(X,Y);
    mse_val = norm(D)/sqrt(num_elements);
end
