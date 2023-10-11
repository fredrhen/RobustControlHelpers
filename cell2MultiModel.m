function [multimodel] = cell2MultiModel(cellArray)

    function multimodel = recursion(index, multimodel, cell_array)
        current_dim = length(index);
        multimodel_dims =  ndims(cell_array);

        if current_dim == multimodel_dims 
            multimodel(:, :, index{:}) = cell_array{index{:}};
        else
            for i=1:size(cell_array, current_dim+1)
            index{current_dim +1 } = i; 
            multimodel = recursion(index, multimodel, cell_array);
            end
        end
        
    end
    
    multimodel = repsys(cellArray{1}, [1 1, size(cellArray)]);    
    multimodel = recursion({}, multimodel, cellArray);
    
end

