function paramstr=GetParamStr(stimulus)
%utility function to generate standardized param strings to use when making stimulus protocols
%note that there is a 255 character limit to the entire emssage string, so
%we will employ some shorthand as needed

paramstr=sprintf('%s',stimulus.type);
fns=fieldnames(stimulus.param);
for fn=1:length(fns)
    if ischar(getfield(stimulus.param, fns{fn}))
        paramstr=[paramstr sprintf(' %s:%s',fns{fn},getfield(stimulus.param, fns{fn}))];
    else
        paramstr=[paramstr sprintf(' %s:%g',fns{fn},getfield(stimulus.param, fns{fn}))];
    end
end

paramstr=strrep(paramstr, 'VarLaser', 'VL');
if length(paramstr)>255
    error('param string is too long. Ask Mike for help')
end

