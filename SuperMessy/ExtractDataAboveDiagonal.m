function data = ExtractDataAboveDiagonal(M)
    % carolina: I made this function, since the one used by the authors was not provided)
    data = [];
    M= triu(M,1);
    N= size(M,1)-1;
    for i = 1:1:N
        data =  cat(2, data, M(i, i+1:end));
    end
end