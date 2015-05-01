
function [level] = otsuthresh(Im_input, pixelid)

% modified otsuthresholding using graythresholding

if ~isempty(pixelid)
    im = Im_input(pixelid);
else
    im = Im_input(:);
end
minval= min(im);
maxval= max(im);
im2 = (im - minval) / (maxval - minval);
level = (minval + (maxval - minval) * graythresh(im2));

end