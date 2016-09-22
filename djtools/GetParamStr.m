function paramstr=GetParamStr(stimulus)
%utility function to generate standardized param strings to use when making stimulus protocols
paramstr=sprintf('%s',stimulus.type);
fns=fieldnames(stimulus.param);
for fn=1:length(fns)
    if ischar(getfield(stimulus.param, fns{fn}))
        paramstr=[paramstr sprintf(' %s:%s',fns{fn},getfield(stimulus.param, fns{fn}))];
    else
        paramstr=[paramstr sprintf(' %s:%g',fns{fn},getfield(stimulus.param, fns{fn}))];
    end
end

