% Zhan Ma @ Thomson Intern
% Yuv ReadIn
% Orignated at Sept. 16. 2008
% Updated at 02/17/2009

%supported format
%yv12: 4:2:0 format
%yv16: 4:2:2 format
%yv24: 4:4:4 format
function [Ymtx, Umtx, Vmtx] = yuv_read(video_src, bit_depth, width, height, frm_offset, pixel_idf)

% video_src: location of video source
% bit_depth: bits per pixel
% width:     picture width
% height:    picture height
% frm_offset: frame offset from the start point, 0 as 1st frame without
% offset
% pixel_idf: could be 'yv12' for 4:2:0, 'yv16' for 4:2:2 and 'yv24' for
% 4:4:4

switch pixel_idf
    case 'yv12'
        pix_idf = 0;
    case 'yv16'
        pix_idf = 1;
    case 'yv24'
        pix_idf = 2;
    otherwise
        sprintf('%c is a wrong pixel sampling patter.\n',pix_idf)
        return
end

%currently, only support 8~16 bit per sample
if ( (bit_depth < 8) || (bit_depth > 16))
    sprintf('Currently, only supporting 8~16 bits per sample, %d is invalid.\n', bit_depth)
    return
end

%currnetly, both width and height should be the factor of 8.
if (( (mod(width,8))) ||(mod(height,8) ))
    sprintf('Invalid width or height of input video, width=%d, height=%d.\n',width, height)
    return
end

%Initilization
if (bit_depth == 8)
   open_mode = 'uint8';
   scal_factor = 1;
else 
   open_mode = 'uint16';
   scal_factor = 2;
end

byte_len   = 1;

if (pix_idf == 0)
    img_size_in_byte = 1.5*width*height*scal_factor;
    Y_size = height*width;
    U_size = Y_size/4;
    V_size = U_size;
    %byte_offset = byte_len*frm_offset*img_size_in_byte;   
elseif (pix_idf == 1)
    img_size_in_byte = 2*width*height*scal_factor;
    Y_size = height*width;
    U_size = Y_size/2;
    V_size = U_size;
else
    img_size_in_byte = 3*width*height*scal_factor;
    Y_size = height*width;
    U_size = Y_size;
    V_size = Y_size;
end

byte_offset= byte_len*frm_offset*img_size_in_byte;
%bit_offset = 1.5*byte_len*frm_offset*width*height*scal_factor;

%Y_size = height*width;
%U_size = Y_size/4;
%V_size = Y_size/4;
%yuv420

Y = zeros(Y_size, 1, open_mode);
U = zeros(U_size, 1, open_mode);
V = zeros(V_size, 1, open_mode);

%input
yuv_fid = fopen(video_src, 'r');
fseek(yuv_fid, 0, 'eof');
video_src_len = ftell(yuv_fid);

%total frame number
tot_frm_num = video_src_len/img_size_in_byte;

if (frm_offset > tot_frm_num)
    sprintf('Invalid Frame Offset %d > Total Frame %d.\n', frm_offset, tot_frm_num);
    return
end
% if (bit_offset > bitstream_len)
%     sprintf('Invalid Frame offset')
%     return
% end

fseek(yuv_fid, byte_offset, 'bof'); % return back pointer to origin
Y       = fread(yuv_fid, Y_size, open_mode);

%byte_offset = byte_offset + byte_len*width*height*scal_factor;
%fseek(yuv_fid, byte_offset, 'bof');

U       = fread(yuv_fid, U_size, open_mode);

%byte_offset = byte_offset + byte_len*(width/2)*(height/2)*scal_factor;
%fseek(yuv_fid, bit_offset, 'bof');

V       = fread(yuv_fid, V_size, open_mode);

fseek(yuv_fid, 0, 'bof'); % return back pointer to the origin

%regulation
if (pix_idf == 0)
    Ymtx = zeros(height, width, open_mode);
    Umtx = zeros(height/2,width/2, open_mode);
    Vmtx = zeros(height/2,width/2, open_mode);
    UV_height = height/2;
    UV_width = width/2;
elseif (pix_idf == 1)
    Ymtx = zeros(height, width, open_mode);
    Umtx = zeros(height, width/2, open_mode);
    Vmtx = zeros(height, width/2, open_mode);
    UV_height = height;
    UV_width = width/2;
else
    Ymtx = zeros(height, width, open_mode);
    Umtx = zeros(height, width, open_mode);
    Vmtx = zeros(height, width, open_mode);
    UV_height = height;
    UV_width = width;
end

for j=1:1:height
    for i=1:1:width
        Ymtx(j,i) = Y(i + width*(j-1));
    end
end
for j=1:1:UV_height
    for i=1:1:UV_width
        Umtx(j,i) = U(i + (UV_width)*(j-1));
    end
end

for j=1:1:UV_height
    for i=1:1:UV_width    
        Vmtx(j,i) = V(i + (UV_width)*(j-1));
    end
end

fclose(yuv_fid);

