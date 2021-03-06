%EPS2PDF  Convert an eps file to pdf format using ghostscript
%
% Examples:
%   eps2pdf source dest
%   eps2pdf(source, dest, crop)
%   eps2pdf(source, dest, crop, append)
%   eps2pdf(source, dest, crop, append, gray)
%   eps2pdf(source, dest, crop, append, gray, quality)
%
% This function converts an eps file to pdf format. The output can be
% optionally cropped and also converted to grayscale. If the output pdf
% file already exists then the eps file can optionally be appended as a new
% page on the end of the eps file. The level of bitmap compression can also
% optionally be set.
%
% This function requires that you have ghostscript installed on your
% system. Ghostscript can be downloaded from: http://www.ghostscript.com
%
%IN:
%   source - filename of the source eps file to convert. The filename is
%            assumed to already have the extension ".eps".
%   dest - filename of the destination pdf file. The filename is assumed to
%          already have the extension ".pdf".
%   crop - boolean indicating whether to crop the borders off the pdf.
%          Default: true.
%   append - boolean indicating whether the eps should be appended to the
%            end of the pdf as a new page (if the pdf exists already).
%            Default: false.
%   gray - boolean indicating whether the output pdf should be grayscale or
%          not. Default: false.
%   quality - scalar indicating the level of image bitmap quality to
%             output. A larger value gives a higher quality. quality > 100
%             gives lossless output. Default: ghostscript prepress default.

% Copyright (C) Oliver Woodford 2009-2010

% Suggestion of appending pdf files provided by Matt C at:
% http://www.mathworks.com/matlabcentral/fileexchange/23629

% Thank you to Fabio Viola for pointing out compression artifacts, leading
% to the quality setting.

function eps2pdf(source, dest, crop, append, gray, quality)
% Intialise the options string for ghostscript
options = ['-q -dNOPAUSE -dBATCH -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile="' dest '"'];
% Set crop option
if nargin < 3 || crop
    options = [options ' -dEPSCrop'];
end
% Set the grayscale option
if nargin > 4 && gray
    options = [options ' -sColorConversionStrategy=Gray -dProcessColorModel=/DeviceGray'];
end
% Set the bitmap quality
if nargin > 5 && ~isempty(quality)
    options = [options ' -dAutoFilterColorImages=false -dAutoFilterGrayImages=false'];
    if quality > 100
        options = [options ' -dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode'];
    else
        options = [options ' -dColorImageFilter=/DCTEncode -dGrayImageFilter=/DCTEncode'];
        v = 1 + (quality < 80);
        quality = 1 - quality / 100;
        s = sprintf('<< /QFactor %.2f /Blend 1 /HSample [%d 1 1 %d] /VSample [%d 1 1 %d] >>', quality, v, v, v, v);
        options = sprintf('%s -c ".setpdfwrite << /ColorImageDict %s /GrayImageDict %s >> setdistillerparams"', options, s, s);
    end
end
% Check if the output file exists
if nargin > 3 && append && exist(dest, 'file') == 2
    % File exists - append current figure to the end
    tmp_nam = tempname;
    % Copy the file
    copyfile(dest, tmp_nam);
    % Add the output file names
    options = [options ' -f "' tmp_nam '" "' source '"'];
    try
        % Convert to pdf using ghostscript
        ghostscript(options);
    catch
        % Delete the intermediate file
        delete(tmp_nam);
        rethrow(lasterror);
    end
    % Delete the intermediate file
    delete(tmp_nam);
else
    % File doesn't exist or should be over-written
    % Add the output file names
    options = [options ' -f "' source '"'];
    % Convert to pdf using ghostscript
    ghostscript(options);
end
return

