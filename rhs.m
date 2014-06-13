function [rep] = rhs(y)
rep = [y(2); y(1)^3-y(1)];
end
