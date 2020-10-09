function imp = implementer(object_cell_array, connection_vector, mode)
    %function implements the connection/disconnection specified by
    %connection_vector and mode (which should be one of 'sb','ub' or 'mb')
    %and returns the resulting objects in a cell array
    
    %first confirm inputs make some sense
    if length(object_cell_array) >1 && (strcmp(mode,'sb') + strcmp(mode,'ub')) == 1
        disp('Error, mode and number of objects passed to implementer do not match')
        return
    end
    [~,ncs] = size(connection_vector);
    if ncs>1
        disp('Error, too many connections passed to implementer')
        return
    end
    
    % branch into various modes
    switch mode
        
        % self-bonding mode
        case 'sb'
            object = object_cell_array{1};
            imp = self_connect(object, connection_vector);
            
        % disconnecting mode    
        case 'ub'
            object = object_cell_array{1};
            imp = disconnect(object, connection_vector);
            
        % object-object bonding mode
        case 'mb'
            
            object1 = object_cell_array{1};
            object2 = object_cell_array{2};
            ele1 = connection_vector{1};
            ele2 = connection_vector{2};
            if ele1(1,1) == 1
            imp = connect_objects(object1, object2, ele1(1,2),ele2(1,2));
            end
            if ele2(1,1) == 1
            imp = connect_objects(object1, object2, ele2(1,2), ele1(1,2));
            end
    end %switch

end %imp