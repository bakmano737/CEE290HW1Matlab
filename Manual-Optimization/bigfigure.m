function varargout = bigfigure(nRows,nCols,figIndex,varargin)

if nargin==0
    nRows = 1;
    nCols = 1;
    figIndex = 1;
end


screenRect = get(0,'screensize');
f = 0.4;
hw = floor(f*screenRect(3));
hh = floor(f*screenRect(4));
mw = floor(0.5*screenRect(3));
mh = floor(0.5*screenRect(4));
rect = [mw-hw,mh-hh,2*hw,2*hh];
clear f hw hh mw mh


border = [86,9,0,9];

existingFigures = get(0,'children');
figureNumber = 1;
while true
    if ismember(figureNumber,existingFigures)
        figureNumber = figureNumber+1;
    else
        break
    end
end

authorizedOptions ={'rect','border','figureNumber'};


for k = 1:2:length(varargin(:))
    if ~strcmp(varargin{k}, authorizedOptions)
        s=dbstack;
        error(['Unauthorized parameter name ' 39 varargin{k} 39 ' in ' 10,...
            'parameter/value passed to ' 39 s(end-1).name 39 '.']);
    end
    eval([varargin{k},'=varargin{',num2str(k+1),'};'])
end

if ~isempty(varargin)
    clear k
end

curRectWidth = floor(rect(3)/nCols);
curRectHeight = floor(rect(4)/nRows);


colNumber = mod(figIndex-1,nCols)+1;
rowNumber = 1+(figIndex-colNumber)/nCols;

if ~(1<=colNumber && colNumber<=nCols)||...
   ~(1<=rowNumber && rowNumber<=nRows)

    error('Figure number out of bounds.')
    
end

curRectLeft = rect(1)+((colNumber-1)*curRectWidth)+1;
curRectBottom = rect(2)+(nRows-rowNumber)*curRectHeight+1;


h=figure(figureNumber);
set(h,'position',[curRectLeft + border(4),...
                     curRectBottom + border(3),...
                     curRectWidth - border(2) - border(4),...
                     curRectHeight - border(1) - border(3)]);



if nargout==1
    varargout={h};
end
                 
                 