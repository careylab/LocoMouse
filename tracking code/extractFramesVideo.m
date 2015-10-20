function varargout = extractFramesVideo(varargin)
% EXTRACTFRAMESVIDEO Writes the frames of a given video file into images.
%
% Usage:
%
% extractFramesVideo(file, output_folder, output_format, options);
%
% Input:
% file: video file path.
% output_folder: output folder path.
% output_format: sprtinf template format for writing the images, including
% image format (e.g. image_%06d.png for image_000001.png etc.).
% options:
%
% 'Background': path to or background image.
% 'Flip': logical value to check for mouse running from left or right.
% 'Undistort': Undistortion map for the bottom half of the image.
%
% example:
%
% extractFramesVideo('C:\example.avi', 'C:\', 'images_%06d.png', 'Background',B,'Undistort',IDX,'Flip',true);
%
% This would read example.avi and save it in C:\ as images_000001.png, etc.
% but first removing the background stored in image B, undistorting the
% image acoording to IDX and checking if the image needed to be flipped
% vertically.

file = varargin{1};
output_folder = varargin{2};
output_format = varargin{3};
st = 3;
Bkg = [];
flip = false;
rescale = false;
crop = false;
rep_option = true(1,5);
vid = VideoReader(file);
mirror_line = vid.Height;
IDX = 1:(mirror_line*vid.Width);
crop_window = [1 1 vid.Height vid.Width];
scale = 1;

while st < length(varargin)
    opt_str = varargin{st+1};
    
    if ~ischar(opt_str)
        error('Options must be specified by strings. Allowed options are ''Flip'', ''Background'', ''Undistort'' and ''Rescale''.');
    end
    
    switch lower(varargin{st+1})
        case 'background'
            if rep_option(1)
                switch class(varargin{st+2})
                    case 'uint8'
                        Bkg = varargin{st+2};
                    case 'char'
                        if strcmpi(varargin{st+2},'median')
                            Bkg = read(vid);
                            Bkg = median(Bkg(:,:,1,:),4);
                        else
                            Bkg = imread(varargin{st+2});
                        end
                    otherwise
                        error('Option ''Background'' must be an image or a path to an image.');
                end
                st = st+2;
                rep_option(1) = false;
            else
                error('Background option defined more than once.');
            end
        case 'flip'
            if rep_option(2)
                if islogical(varargin{st+2})
                    flip = varargin{st+2};
                else
                    error('Option ''Flip'' must contain a boolean');
                end
                st = st+2;
                rep_option(2) = false;
            else
                error('Flip option defined more than once');
            end
        case 'undistort'
            if rep_option(3)
                IDX = varargin{st+2}{1};
                mirror_line = varargin{st+2}{2};
            else
                error('Undistort option defined more than once');
            end
            st = st+2;
            rep_option(3) = false;
        case 'rescale'
            if rep_option(4)
                scale = varargin{st+2};
            else
                error('Rescale option defined more than once!');
            end
            rescale = true;
            st = st+2;
            rep_option(4) = false;
        case 'crop'
            if rep_option(5)
                crop_window = varargin{st+2};
            else
                error('Crop option defined more than once!');
            end
            
            crop = true;
            st = st+2;
            rep_option(5) = false;
            
        otherwise
            error('Unknown option. Allowed options are ''Flip'' and ''Background''.');
    end
end


Nframes = get(vid,'NumberOfFrames');

if flip
    % flip being on does not mean to flip, but to look for a flip!
    flip = checkMouseSide(vid,Bkg,10); 
end

if isempty(Bkg)
    check_bkg = false;
else
    check_bkg = true;
end

repoption3 = rep_option(3);

if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

parfor i = 1:Nframes
    I = read(vid,i);
    %% for led activation
    I2 = I;
    if ndims(I) > 2
        I = rgb2gray(I);
    end
    
    if check_bkg
        I2 = I-Bkg;
    end
    
    if ~repoption3
        if numel(IDX) == numel(I)
            I2(:) = I2(IDX);
        else
            [I_top,I_bottom] = splitImage(I2,mirror_line);
            N = uint8(zeros(size(I_bottom)));
            N(:) = I_bottom(IDX(:));
            I2 = [I_top;N];
        end
    end
    
    if crop
        cw = crop_window;
        I_res = I2(:,cw(1)+cw(2)+1:end);
        I2 = I2(:,cw(1)-1 + (1:cw(2)));
    else
        I_res = [];
    end
    
    if flip
        I2 = I2(:,end:-1:1);
    end
    
    if rescale
        I2 = imresize(I2,scale);
    end
    
    imwrite(I2,fullfile(output_folder,sprintf(output_format,i)));
    if crop
        imwrite(I_res,fullfile(output_folder,sprintf(['res_image_' output_format],i)));
    end
end

% if needed output flip status:
if nargout > 0
    varargout{1} = flip;
end