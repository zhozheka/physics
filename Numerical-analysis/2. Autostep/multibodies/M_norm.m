function [ N ] = M_norm( Matrix )
    N = sqrt(sum(Matrix(:).^2))
end

