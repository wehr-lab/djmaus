
function Message(text,type)
global SP
h = SP.Messageh;

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
        
        if size(prev_text,1)>=6
            new_text=sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s',prev_text(end-5,:),prev_text(end-4,:),prev_text(end-3,:),prev_text(end-2,:),prev_text(end-1,:),prev_text(end,:), text);
        elseif size(prev_text,1)==5
            new_text=sprintf('%s\n%s\n%s\n%s\n%s\n%s',prev_text(end-4,:),prev_text(end-3,:),prev_text(end-2,:),prev_text(end-1,:),prev_text(end,:), text);
        elseif size(prev_text,1)==4
            new_text=sprintf('%s\n%s\n%s\n%s\n%s',prev_text(end-3,:),prev_text(end-2,:),prev_text(end-1,:),prev_text(end,:), text);
        elseif size(prev_text,1)==3
            new_text=sprintf('%s\n%s\n%s\n%s',prev_text(end-2,:),prev_text(end-1,:),prev_text(end,:), text);
        elseif size(prev_text,1)==2
            new_text=sprintf('%s\n%s\n%s',prev_text(end-1,:),prev_text(end,:), text);
        elseif size(prev_text,1)==1
            new_text=sprintf('%s\n%s',prev_text, text);
        elseif isempty(prev_text)
            new_text= text;
        else
            error('??????')
        end
        
        try
            set(h,'string',new_text,'backgroundcolor',[1 1 1], 'foregroundcolor', [0 0 0]);
        catch
            set(h,'string',text,'backgroundcolor',[1 1 1], 'foregroundcolor', [0 0 0]);
        end
    otherwise %no type requested
        set(h,'string',text,'backgroundcolor',get(gcf,'color'));
end
drawnow
