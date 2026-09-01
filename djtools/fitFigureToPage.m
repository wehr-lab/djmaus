function fitFigureToPage(fig, pageW, pageH, margin)
%FITFIGURETOPAGE resizes fig to be exactly pageW x pageH inches (e.g. for
%exportgraphics), rescaling and repositioning its axes/colorbars/legends
%as a group so the actual plotted content -- not the original figure
%canvas, which often has extra unused space around the content -- is
%scaled to fit within the page (minus margin) without stretching, and
%centered. Reproduces print's '-bestfit' behavior for exportgraphics,
%which has no such option itself.
%
%usage: fitFigureToPage(fig, 8.5, 11, 0.25)  %US Letter portrait, .25in margin (matches print's -bestfit)
%       fitFigureToPage(fig, 11, 8.5)        %US Letter landscape, default margin
%
%Note: this repositions axes/subplots, colorbars, and legends (found via
%findall). Textbox/annotation objects with their own independent Position
%(created via `annotation(...)`, not axes text) are NOT covered and would
%need the same scale/offset treatment added if a figure uses them.
%-mike/claude 9.1.26

if nargin<4, margin=0.25; end %matches print's -bestfit margin

fig.Units='inches';
origFigPos=fig.Position; %[left bottom width height] of the figure as originally laid out

objs=[findall(fig, 'type', 'axes'); findall(fig, 'type', 'colorbar'); findall(fig, 'type', 'legend')];
if isempty(objs)
    fig.Position=[origFigPos(1) origFigPos(2) pageW pageH]; %nothing to rescale -- just resize the figure itself
    return
end

%Capture every repositionable child's position in inches BEFORE resizing
%the figure. This matters because Position is commonly stored in
%normalized units (fraction of the parent figure's current size) --
%converting Units to 'inches' re-expresses Position relative to whatever
%size the figure is AT THAT MOMENT. If we resized the figure first, this
%conversion would be relative to the new page size instead of the
%original figure, and our subsequent scale/offset math would be applied
%on top of an already-wrong baseline.
origPos=cell(numel(objs),1);
for i=1:numel(objs)
    objs(i).Units='inches';
    origPos{i}=objs(i).Position; %[left bottom width height], relative to the ORIGINAL figure's coordinate frame
end

%Compute the tight bounding box enclosing all objects' content, rather
%than assuming content spans the whole original figure canvas -- several
%figures in this codebase have substantial unused canvas around their
%actual plotted content, and preserving THAT ratio (instead of the tight
%content) is what was producing pages that were mostly blank with tiny
%content shoved off to one side.
lefts   = cellfun(@(p) p(1), origPos);
bottoms = cellfun(@(p) p(2), origPos);
rights  = cellfun(@(p) p(1)+p(3), origPos);
tops    = cellfun(@(p) p(2)+p(4), origPos);
bboxLeft=min(lefts); bboxBottom=min(bottoms);
contentW=max(rights)-bboxLeft;
contentH=max(tops)-bboxBottom;

availW=pageW-2*margin;
availH=pageH-2*margin;
scale=min(availW/contentW, availH/contentH); %preserves aspect ratio, no stretching

offsetX=(pageW-contentW*scale)/2;
offsetY=(pageH-contentH*scale)/2;

fig.Position=[origFigPos(1) origFigPos(2) pageW pageH]; %now resize the figure itself to the full page

for i=1:numel(objs)
    p=origPos{i};
    %reposition relative to the content bounding box's own origin, not
    %the whole original figure's origin
    objs(i).Position=[offsetX+(p(1)-bboxLeft)*scale, offsetY+(p(2)-bboxBottom)*scale, p(3)*scale, p(4)*scale];
end
end
