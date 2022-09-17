function struct_out = StructToDouble(struct_in)
%change all variables in struct_in (which is a cell of structs) to single/double
%precision

struct_out = struct_in;

for i = 1:length(struct_out)
    fn = fieldnames(struct_out{i});
    for j = 1:length(fn)
        if ~(strcmp(fn{j},'distance_spectrum') || strcmp(fn{j},'joint_type'))
            struct_out{i}.(fn{j}) = single(struct_out{i}.(fn{j}));
        end
    end
end

end

