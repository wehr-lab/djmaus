% example of Report Generator

%create 3 figures
close all
figure
hold on
for i=1:20
    plot(randn(1, 100)+i)
end
xlabel('stuff')
set(gcf, 'pos', [   680    84   551   894])
figure
plot(randn(5))
xlabel('stuff')
figure
hist(randn(5))
xlabel('stuff')
title('distribution')

% Import the base classes
import mlreportgen.report.*
import mlreportgen.dom.*

% Create a report object
% rpt = Report('example', 'html');
 rpt = Report('example', 'pdf');

% Add a title page
% To customize additional title page properties, see mlreportgen.report.TitlePage.
tp = TitlePage;
tp.Title = 'MouseID_Experimentat';
tp.Subtitle = 'various stuff';
tp.Author = 'Mike';
add(rpt,tp);

% Add default table of contents
add(rpt,TableOfContents);
chap=Chapter;
    chap.Title = 'this is a chapter';
para=Paragraph ('here is some text');
add(chap, para);
add(rpt, chap);

for i=1:5
    % Add a section
    ch(i) = Section;
    ch(i).Title = sprintf('cell %d', i);
    para = Paragraph(sprintf('this is cell number %d', i)) ;
    add(ch(i),para)
    
    
    % Add a figure.
    % For more information on figures, see mlreportgen.report.Figure.
    % https://www.mathworks.com/help/rptgen/ug/mlreportgen.report.figure-class.html
    
    for j=1:3
        fig(j)=Figure(figure(j));
        figi(j) = Image(getSnapshotImage(fig(j),rpt));
    end
    
    table=Table(2);
    table.Border = 'single';
    table.ColSep = 'single';
    table.RowSep = 'single';
    row = TableRow;
    te = TableEntry(figi(1));
%    te = TableEntry('1');
    te.RowSpan = 2;
    te.ColSpan = 1;
    append(row, te);
    append(row, TableEntry(figi(2)));
%     append(row, TableEntry('2'));
    append(table,row);
    row = TableRow;
   append(row, TableEntry(figi(3)));
%     append(row, TableEntry('3'));
    append(table,row);
    
    
    add(ch(i),table);
    add(rpt,ch(i))
    
end

% Close and run the report.

close(rpt)
rptview(rpt)



