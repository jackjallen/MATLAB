function Y = str2strarray(X)

Y = [strrep(fieldnames(X),'.','_'),  cellfun(@(x) iff(isnumeric(x),num2cell(x),x),struct2cell(X),'un',0)].';

Y = struct(Y{:});
