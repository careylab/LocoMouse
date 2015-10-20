function [I_top, I_bottom] = splitImage(I, mirror_line)
% SPLITIMAGE    Splits the image into upper and lower half according to the
% level set by mirror_line.
%
% Input:
% 
% I: the NxM image.
% mirror_line: an integer between 0 and N. 0 returns I in I_bottom, N
% returns I in I_top.

N = size(I,1);

if mirror_line < 0 || mirror_line > N
    error('mirror_line must be between 0 and N');
end

mirror_line = round(mirror_line);

I_top = I(1:mirror_line,:);
I_bottom = I(mirror_line+1:end,:);