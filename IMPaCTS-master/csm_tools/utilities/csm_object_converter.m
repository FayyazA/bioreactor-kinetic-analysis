classdef csm_object_converter
%CSM_OBJECT_CONVERTER - Convert an object into a struct, and vice versa;
%
% Usage:
%
%	csm_struct = csm_object_converter.toStruct( csm_object );
%	csm_object = csm_object_converter.toObject( csm_struct );
%
% Methods:
%
%	csm_object_converter.toStruct( csm_object ) : Convert an object to a struct.
%	csm_object_converter.toObject( csm_struct ) : Convert an exported struct into an object.
%
% Description:
%	
%	This class contains methods for exporting Matlab objects to structs,
%	and vice versa.
%
%	This class is necessary because if you share data to someone who does
%	not have the toolbox, or a class definition has been updated, the old
%	object will not import to their workspace.
%
%	Does not support Cell arrays containing objects. ie it will not
%	recursively loop into a cell array. Convert cell to struct before
%	using.
%
% Copyright (C) Division of Computational and Systems Medicine, Imperial College London - 2014

% Author - Gordon Haggart 2016


    methods( Static )

        % Recursively create structs from objects.
        % If a property of an object is another object, it loops into that
        % object as well.
        %
        % Uses evals because of dynamic property names
        % Creates csm_class_info for mapping variable names to object
        % types.
        function [output] = toStruct( object, no_warning )

            if nargin < 2 || false( no_warning )

                warning( 'This method does not support Cell arrays containing objects! See description for more information.');

            end

            output = struct;

            % If its a normal datatype (Single), just return the value
            if csm_object_converter.isSingleDatatype( object )

                output = object;

                return;

            end

            % Get the properties

            if isa( object, 'struct' )

                property_list = fieldnames( object );

            elseif isa( object, 'containers.Map' )

                property_list = keys( object );

            else

                property_list = properties( object );
                
                output.csm_class_type = class( object );
            
                output.csm_subclass_info = containers.Map;
                
            end

            for i = 1 : length( property_list )

                property_accessor = csm_object_converter.buildPropertyAccessor( object, property_list, i  );

                % If its a standard data type, just add it to the output struct
                if csm_object_converter.isSingleDatatype( property_accessor )

                    output.( property_list{ i } ) = property_accessor;


                % If its a struct, create a new struct, loop through and add 1 by 1.
                elseif isa( property_accessor,  'struct' )
                    
                    new_struct = struct;

                    old_struct = property_accessor;

                    for p = 1 : length( old_struct )
                    
                        new_struct = csm_object_converter.toStruct( old_struct( p ), true  );
                        
                    end 
                    
                    output.( property_list{ i } ) = new_struct;
                                
                    
                % If its a map container, create a new map container, loop through and add 1 by 1
                elseif isa( property_accessor,  'containers.Map' )
                    
                    new_map = containers.Map;

                    old_map = property_accessor;
                    
                    map_keys =  keys( old_map );
                    
                    for p = 1 : length( map_keys )

                        map_value = old_map( map_keys{ p } );

                        new_map( map_keys{ p } ) = csm_object_converter.toStruct( map_value, true );
                        
                    end
                    
                    output.( property_list{ i } ) = new_map;
                
                else

                    if isfield(output,'csm_subclass_info')
                    
                        % Save the object type in the csm_class_info cell array.
                        output.csm_subclass_info( property_list{ i } ) = class( property_accessor );
                    end    

                    % Recall the toStruct method in a nested way
                    output.( property_list{ i } ) = csm_object_converter.toStruct( property_accessor, true  );

                end

            end

        end

        % Check if this is a CSM Toolbox Class
        % Simply examines the object, and whether it is a standard data type
        % Uses evals because of the dynamic object property names.
        %
        % Sample code used:
        %	if strcmp( class( object.property ), 'int8' )
        %		it_is = false;
        %	end
        %
        function [it_is] = isSingleDatatype( variable )

            it_is = false;

            if isa( variable, 'int8' )

                it_is = true;

            end

            if isa( variable, 'int16' )

                it_is = true;

            end

            if isa( variable, 'uint16' )

                it_is = true;

            end

            if isa(  variable, 'int32' )

                it_is = true;

            end

            if isa( variable, 'uint32' )

                it_is = true;

            end

            if isa( variable, 'int64' )

                it_is = true;

            end

            if isa( variable, 'uint64' )

                it_is = true;

            end

            if isa( variable, 'single' )

                it_is = true;

            end

            if isa( variable, 'double' )

                it_is = true;

            end

            if isa( variable, 'logical' )

                it_is = true;

            end

            if isa( variable, 'char' )

                it_is = true;

            end

            if isa( variable, 'cell' )

                it_is = true;

            end
            
            if isa( variable, 'dataset' )
                
                it_is = true;
                
            end  
            
            if isa( variable, 'table' )
                
                it_is = true;
                
            end  
        
        end
        
        function [it_is] = inIgnoreList( structure )
            
            it_is = false;
            
            if isstruct(structure) && isfield(structure,'csm_class_type')
                
               class_name = structure.csm_class_type;
               
               if strcmp(class_name,'inputParser')
                   
                   it_is = true;
                   
               end    
                
            end    
            
        end    

        % Checks if its in the map container, and then whether its a csm class
        function [it_is] = isCsmObject( csm_subclass_info, property_name )

            it_is = false;
            
            if strcmp( property_name, 'inputparser' )
                
                return;
                
            end    

            if isKey( csm_subclass_info, property_name )

                % Check the field exists and is a class
                if (isfield( csm_subclass_info( property_name ), 'class' )) || (exist('csm_audit_info','class') == 8)

                    it_is = true;

                else

                    warning( strcat( 'Property is in csm_subclass_info Map, but is not a class. Check exporter for ', property_name ));

                end

            end

        end
        
        
        function [property_accessor] = buildPropertyAccessor( object, property_list, i )

            if isa( object, 'containers.Map' )

                property_accessor = object( property_list{ i } );
                
            else

                property_accessor = object.( property_list{ i });

            end


        end
        
        function [output] = toObject( structure )
                        
            % If its a normal datatype (Single), just return the value
            if csm_object_converter.isSingleDatatype( structure ) || csm_object_converter.inIgnoreList( structure )

                output = structure;

                return;

            end
            
            % If it's a struct, loop over it like a struct
            if isstruct( structure )
            
                fields = fieldnames( structure );
                
                is_object = false;

                if isfield(structure,'csm_subclass_info')
                    
                    if strcmp(structure.csm_class_type,'csm_import_raw_nmr')
                        
                        eval( strcat( 'output = ', structure.csm_class_type ,'("build_empty_object",true);'));
                    
                    else
                        
                        eval( strcat( 'output = ', structure.csm_class_type ,'();'));
                   
                    end
                    
                    is_object = true;
                
                else 
                    
                    output = struct;
                    
                end    
                    
                for i = 1 : numel( fields )
                    
                    if ~strcmp(fields{i},'csm_subclass_info') && ~strcmp(fields{i},'csm_class_type')

                        temp_output = csm_object_converter.tObject( structure.(fields{i}));
                    
                        disp(fields{i});
                        
                        if strcmp(fields{i},'Quantile')
                            break_here = 1;
                        end    
                         
                        if is_object
                            
                            output = setProperty(output,fields{i},temp_output);
                            
                        else
                           
                            output.(fields{i}) = temp_output;
                            
                        end    
                    end
                end

            elseif isa( structure,  'containers.Map' )   
                
                output = containers.Map;

                map_keys =  keys( structure );

                for p = 1 : length( map_keys )

                    map_value = structure( map_keys{ p } );

                    output( map_keys{ p } ) = csm_object_converter.tObject( map_value );

                end
                
            end    
        end    

    end

end

