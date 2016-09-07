
function djMessage(text,type)
global SP
h = SP.Messageh;

if isempty(text)
    return %don't bother displaying nothing
end
%useful for debugging:
fprintf('\n%s', text)

if nargin < 2 type = ''; end
switch type
    case 'clear'
        set(h,'string','','backgroundcolor',get(findobj('tag',module),'color'));
    case 'error'
        set(h,'string',text,'backgroundcolor',[1 0 0]);
        
    case 'blue'
        set(h,'string',text,'backgroundcolor',[0 0 1], 'foregroundcolor', [1 1 1]);
    case 'green'
        set(h,'string',text,'backgroundcolor',[0 1 0], 'foregroundcolor', [0 0 0]);
    case 'normal'
        set(h,'string',text,'backgroundcolor',get(findobj('tag',module),'color'), 'foregroundcolor', [0 0 0]);
    case 'append' %mw 06.26.06
        prev_text=get(h, 'string');
        if iscell(prev_text)
            new_text={prev_text{:}, text};
        else
            new_text={prev_text, text};
        end
        wrapped_text=textwrap(SP.Messageh, new_text, 40);
        
        try
            set(h,'string',wrapped_text,'backgroundcolor',[1 1 1], 'foregroundcolor', [0 0 0]);
        catch
            set(h,'string',text,'backgroundcolor',[1 1 1], 'foregroundcolor', [0 0 0]);
        end
    otherwise %no type requested
        set(h,'string',text,'backgroundcolor',get(gcf,'color'));
end
drawnow
