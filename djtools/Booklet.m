function out = Booklet(action, pdffile, varargin)
%Booklet - build a multi-page pdf "booklet" of figures as compressed JPEG pages: small files that
%browse quickly over the network. Pages are exported one at a time to a local temp folder and
%written into one pdf at the end, by this function itself (no Report Generator, no external tools).
%
%   Booklet('start',  pdffile)              %begin (clears any leftover pages for this pdf)
%   Booklet('add',    pdffile, fig)         %add a figure as the next page
%   Booklet('add',    pdffile, fig, 'Resolution', 150, 'Quality', 85)
%   n = Booklet('finish', pdffile)          %write pdffile (replacing it), delete the temp pages
%
%Each page is exported with exportgraphics to a lossless PNG (so it is always a raster, whether
%the figure is visible or not), then recompressed to JPEG with imwrite 'Quality'. The page size in
%the pdf is the figure's size in inches (pixels / Resolution), so zooming in a pdf viewer works as
%for the old booklets. Typical page: 150-400 KB at 150 dpi, vs ~20 MB for a vector page with rasters.
%
%name-value options for 'add':
%   'Resolution'  150   dpi of the exported page
%   'Quality'     85    JPEG quality (1-100)
%-mike/claude 10.8.26

action=lower(char(action)); pdffile=char(pdffile);
tdir=pagedir(pdffile);
switch action
    case 'start'
        if isfolder(tdir), rmdir(tdir, 's'); end
        mkdir(tdir); out=tdir;
    case 'add'
        p=inputParser; addRequired(p, 'fig'); addParameter(p, 'Resolution', 150, @isscalar); addParameter(p, 'Quality', 85, @isscalar);
        parse(p, varargin{:}); P=p.Results;
        if ~isfolder(tdir), mkdir(tdir); end
        k=numel(dir(fullfile(tdir, 'p*.jpg')))+1;
        png=fullfile(tdir, sprintf('tmp%05d.png', k)); jpg=fullfile(tdir, sprintf('p%05d.jpg', k));
        exportgraphics(P.fig, png, 'Resolution', P.Resolution);
        im=imread(png); if size(im,3)==1, im=repmat(im, 1, 1, 3); end
        imwrite(im, jpg, 'Quality', P.Quality); delete(png);
        fid=fopen(fullfile(tdir, sprintf('p%05d.dpi', k)), 'w'); fprintf(fid, '%g', P.Resolution); fclose(fid);
        out=k;
    case 'finish'
        d=dir(fullfile(tdir, 'p*.jpg')); [~, o]=sort({d.name}); d=d(o);
        out=numel(d);
        if out==0, warning('Booklet:empty', 'no pages for %s', pdffile); return, end
        dpi=zeros(out,1);
        for k=1:out
            f=fullfile(tdir, strrep(d(k).name, '.jpg', '.dpi'));
            if isfile(f), dpi(k)=str2double(fileread(f)); else, dpi(k)=150; end
        end
        writepdf(pdffile, fullfile(tdir, {d.name}), dpi);
        rmdir(tdir, 's');
    otherwise
        error('Booklet: unknown action ''%s'' (start, add, finish)', action);
end
end

%% ---------------------------------------------------------------------
function tdir=pagedir(pdffile)
%local temp folder for one booklet's pages, keyed by the pdf path
key=regexprep(pdffile, '[^A-Za-z0-9]', '_'); if numel(key)>120, key=key(end-119:end); end
tdir=fullfile(tempdir, 'Booklet', key);
end

function writepdf(pdffile, jpgs, dpi)
%minimal PDF 1.4: catalog, page tree, and per page a page object, a content stream that draws the
%image full-page, and the JPEG itself as a DCTDecode image XObject (embedded unchanged)
n=numel(jpgs); nobj=2+3*n;
buf=cell(nobj,1);
kids=sprintf('%d 0 R ', 3+3*(0:n-1));
buf{1}=sprintf('1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n');
buf{2}=sprintf('2 0 obj\n<< /Type /Pages /Kids [%s] /Count %d >>\nendobj\n', strtrim(kids), n);
for k=1:n
    fid=fopen(jpgs{k}, 'r'); J=fread(fid, Inf, '*uint8')'; fclose(fid);
    [w, h, nc]=jpegsize(J);
    cs='/DeviceRGB'; if nc==1, cs='/DeviceGray'; elseif nc==4, cs='/DeviceCMYK'; end
    pw=w*72/dpi(k); ph=h*72/dpi(k);
    io=3+3*(k-1); co=io+1; mo=io+2;
    buf{io}=sprintf('%d 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 %.2f %.2f] /Resources << /XObject << /Im0 %d 0 R >> >> /Contents %d 0 R >>\nendobj\n', ...
        io, pw, ph, mo, co);
    c=sprintf('q %.2f 0 0 %.2f 0 0 cm /Im0 Do Q', pw, ph);
    buf{co}=sprintf('%d 0 obj\n<< /Length %d >>\nstream\n%s\nendstream\nendobj\n', co, numel(c), c);
    hdr=sprintf('%d 0 obj\n<< /Type /XObject /Subtype /Image /Width %d /Height %d /ColorSpace %s /BitsPerComponent 8 /Filter /DCTDecode /Length %d >>\nstream\n', ...
        mo, w, h, cs, numel(J));
    buf{mo}=[uint8(hdr) J uint8(sprintf('\nendstream\nendobj\n'))];
end
head=uint8(sprintf('%%PDF-1.4\n%%%s\n', char([226 227 207 211])));
off=zeros(nobj,1); pos=numel(head);
for k=1:nobj
    if ischar(buf{k}), buf{k}=uint8(buf{k}); end
    off(k)=pos; pos=pos+numel(buf{k});
end
xref=[sprintf('xref\n0 %d\n0000000000 65535 f \n', nobj+1) sprintf('%010d 00000 n \n', off) ...
      sprintf('trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF\n', nobj+1, pos)];
tmp=[pdffile '.writing'];
fid=fopen(tmp, 'w'); if fid<0, error('Booklet: cannot write %s', tmp); end
fwrite(fid, head); for k=1:nobj, fwrite(fid, buf{k}); end; fwrite(fid, uint8(xref)); fclose(fid);
movefile(tmp, pdffile, 'f');   %write to a temp name first so a crash never leaves a half-written booklet
end

function [w, h, nc]=jpegsize(J)
%width, height, components from the JPEG start-of-frame marker (SOF0-SOF15, except DHT/JPG/DAC)
i=3;
while i<numel(J)
    if J(i)~=255, i=i+1; continue, end
    m=J(i+1);
    if m>=192 && m<=207 && ~ismember(m, [196 200 204])
        h=double(J(i+5))*256+double(J(i+6)); w=double(J(i+7))*256+double(J(i+8)); nc=double(J(i+9)); return
    end
    if m==216 || m==1 || (m>=208 && m<=215), i=i+2; continue, end
    i=i+2+double(J(i+2))*256+double(J(i+3));
end
error('Booklet: no JPEG frame header found');
end
