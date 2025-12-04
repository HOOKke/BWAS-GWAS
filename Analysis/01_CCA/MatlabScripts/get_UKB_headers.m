function colnames=get_UKB_headers(data)

colnames        = data.Properties.VariableDescriptions;
colnames        = regexp(colnames, '''(.+)''', 'tokens', 'once');
empty           = cellfun(@isempty, colnames);
colnames(empty) = data.Properties.VariableNames(empty);
colnames        = vertcat(colnames{:});
