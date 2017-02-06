function [Nlam] = CheckInput(wing_data,Nlam,varargin)

if length(varargin) == 1
    morph = varargin{1};
    morphflag = 1;
else
    morphflag = 0;
end

% -- Check if Spar are defined between skin
skinLE = wing_data(:,9);
skinTE = wing_data(:,10);
Ncol  = size(wing_data,2);
for iCol = 12:Ncol
    if ~isempty(find(wing_data(:,iCol)<skinLE,1)) 
        error(['WingBox definition incorrect. WingData Colum: ',num2str(iCol),' Row: ',num2str(find(wing_data(:,iCol)<skinLE,1))])
    end
    if ~isempty(find(wing_data(:,iCol)>skinTE,1))
        error(['WingBox definition incorrect. WingData Colum: ',num2str(iCol),' Row: ',num2str(find(wing_data(:,iCol)>skinTE,1))])
    end
end


if morphflag == 1 && morph.inp.span.flag == 1
    % Check whether the inboard, double, and outboard section each have
    % only one laminate otherwise throw warning
    if Nlam.Span(morph.inp.span.fixedsec)>1 || Nlam.Span(morph.inp.span.doublesec)>1 || Nlam.Span(morph.inp.span.extsec)>1
        warning(['Please note that in case of span extension, all laminates in the span extension region will be ignored and'...
            ' 1 laminate is used for the inboard section and 1 laminate for the outboard section'])
    elseif Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.fixedsec)))>1 || Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.doublesec)))>1 || Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.extsec)))>1
        warning(['Please note that in case of span extension, all laminates in the span extension region will be ignored and'...
            ' 1 laminate is used for the inboard section and 1 laminate for the outboard section'])
    end
        
    if Nlam.Span(morph.inp.span.extsec)>1
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.extsec-1))+(2:Nlam.Span(morph.inp.span.extsec))) = [];
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.extsec-1))+1) = 1;
    end
    if Nlam.Span(morph.inp.span.doublesec)>1
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.doublesec-1))+(2:Nlam.Span(morph.inp.span.doublesec))) = [];
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.doublesec-1))+1) = 1;
    end
    if Nlam.Span(morph.inp.span.fixedsec)>1
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.fixedsec-1))+(2:Nlam.Span(morph.inp.span.fixedsec))) = [];
        Nlam.Chord(sum(Nlam.Span(1:morph.inp.span.fixedsec-1))+1) = 1;
    end
    Nlam.Span([morph.inp.span.fixedsec,morph.inp.span.doublesec,morph.inp.span.extsec]) = 1;
end
% --- Check if size Nlam is consistent with number of fixed nodes in WingData

if length(Nlam.Span)~=(length(find(wing_data(:,11)))-1)
    error('Length of Nlam.Span must be equal to the number of fixed nodes-1  (Col. 11 in Wingdata)')
end

if sum(Nlam.Span) ~= length(Nlam.Chord)
    error('sum(Nlam.Span) must be equal to the Length of Nlam.Chord')
end




end

