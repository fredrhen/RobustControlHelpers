function [cell_array] = multiModel2Cell(MM)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    function cell_array = recursion(index, multimodel, cell_array)
        current_dim = length(index);
        multimodel_dims =  ndims(multimodel) - 2;

        if current_dim == multimodel_dims 
            cell_array{index{:}} = multimodel(:, :, index{:});
        else
            for i=1:size(multimodel, current_dim+3)
            index{current_dim +1 } = i; 
            cell_array = recursion(index, multimodel, cell_array);
            end
        end
        
    end
    
    multi_size = size(MM);
    multi_size = multi_size(3:end);
    cell_array = cell(multi_size);

    cell_array = recursion({}, MM, cell_array);
    
end

