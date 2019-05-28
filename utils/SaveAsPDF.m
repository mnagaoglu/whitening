function isSuccess = SaveAsPDF(fh, filename)
% isSuccess = SaveAsPDF(fh, filename)
%
% MNA 5/14/19 wrote it.
%

if nargin < 1 || isempty(fh)
    fh = gcf;
end

if nargin < 2 || isempty(filename)
    error('SaveAsPDF:filename is missing or empty.')
end

try
    set(fh,'Units','inches');
    pos = get(fh,'position'); 
    set(gcf,'PaperPositionMode','Auto','PaperUnits','inches',...
        'PaperSize',[pos(3), pos(4)]);
    print(gcf,filename,'-dpdf','-r0')
    isSuccess = true;

catch err
    isSuccess = false;
    err.message
    err.stack.line
    err.stack.name
end