% ===============================================================================================
% NG_Codec Project: Global Motion Based Video Coding
% Video Lab @ New York University, School of Engineering
%
% This Program will ...
% Inputs: (NOW ONLY RETURNS MASK, NEEDS MODIFICATION FOR THE LINES DOWN)
% 1. Current Frame: "Y_Frame"
% 2. Sparse Foreground Image: "E"
% 3. Parameter Configuration: "parameters" a struct containing the following fields
%    - QP: quantization parameter of the frame
%    - Debugging_Enabled: Enable debug flag
%    - Displacement: displacement derived from previous tau, for deriving the bounding box margin
%    - Conn_percent: connectivity area percentage relative to the framesize
%
% Outputs:
% 1. Current SubFrame(s): "Y_SubFrame"
% 2. Bounding Box Information (x, y, width, height): "bounding_box_info"
% 3. Binary Foreground Mask in the Bounding Box: "Mask"
% * in the case of multiple objects, all three outputs are cell arrays.
% ===============================================================================================

function Mask = postprocessing(Y_Frame, E, parameters)

if isfield(parameters,'ByPassFilling'),
    ByPassFilling = parameters.ByPassFilling;
else
    ByPassFilling = 0;
end

if isfield(parameters,'QuadTreeMode'),
    QuadTreeMode = parameters.QuadTreeMode;
else
    QuadTreeMode = 0;
end

if isfield(parameters,'ConfThrLo'),
    fg_thrLo = parameters.ConfThrLo;
else
    fg_thrLo = 6;
end

if isfield(parameters,'ByPassAddMorpFilt'),
    ByPassAddMorpFilt = parameters.ByPassAddMorpFilt;
else
    ByPassAddMorpFilt = 1;
end

% Input frame dimension
[Row, Col] = size(Y_Frame);
% Constants
% morphological_processing_window_size = 0;
connectivity_neighbor_num = 4;
bb_marginLo = 5;
if ~QuadTreeMode,
    stripe_removal_threshold = 8; 	% Remove narrow stripe
%     stripe_removal_threshold = 3;
else
    stripe_removal_threshold = 0;
end
channel = parameters.channel; % To handle U and V with smaller threshold but still adaptive
if channel ~= 1
	fg_thrLo = 5;
end


% Parameters
fg_threshold = max([fg_thrLo, 2^((parameters.QP-4)/6)/2]);
bb_margin = max([bb_marginLo, 1.2*parameters.Displacement]);
IsDebugging = parameters.Debugging_Enabled;
% connectivity_area = round(parameters.Conn_percent * (Row*Col)); % 0.01
connectivity_area = parameters.Conn_area;


% Confidence Thresholding
Ethr = E;
Ethr((abs(E)<=fg_threshold)) = 0;
Ethr((abs(E)>fg_threshold)) = 255;
if IsDebugging
	close all;
	figure;imshow(E,[]);title('TP1 - Original Error Image'); close all;
	figure;imshow(Ethr);title('Binary Sparse Image After Confidence Thresholding');
	print(gcf,'-dpng','Binary_Sparse_Image_After_Confidence_Thresholding.png');
end

% Image Morphological Processing
if ~ByPassAddMorpFilt,
    E3 = imclose(Ethr,strel('disk',2));
%     E3 = imclose(Ethr,strel('square',3));
else
	E3 = Ethr;
end

if IsDebugging
	figure;imshow(E3,[]);title('Binary Sparse Image After Morphological Processing');
	print(gcf,'-dpng','Binary_Sparse_Image_After_Morphological_Processing.png');
end

% Image Connectivity Processing
E4 = bwareaopen(E3, connectivity_area);
if IsDebugging
	figure;imshow(E4,[]);title('Binary Error Image After Connectivity Filtering');
	print(gcf,'-dpng','Binary_Sparse_Image_After_Connectivity_Filtering.png');	
end
% Remove the boundary strips
E4(1:stripe_removal_threshold,:) = 0;
E4(:,1:stripe_removal_threshold) = 0;
E4(end-stripe_removal_threshold+1:end,:) = 0;
E4(:,end-stripe_removal_threshold+1:end) = 0;

E5 = bwlabel(E4,connectivity_neighbor_num);

% Region Extraction and Bounding Box Generation
% coordinate = regionprops(E5,'BoundingBox');
% subMaskImage = regionprops(E5,'Image');
MaskIdxList = regionprops(E5,'PixelIdxList');
% Return if no object
if isempty(MaskIdxList),
    disp('No object found!');
    Mask = [];
    return;
end

if ~ByPassFilling,
    % Hole Filling
    E6 = imfill(im2bw(E5),'hole');
    if IsDebugging
        figure;imshow(E6,[]);title('before 2nd imopen');
    end
    E6 = imerode(E6,strel('square',5));
    
    if IsDebugging
        figure; imshow(E6,[]); title('Binary Error Image After Hole Filling');
        print(gcf,'-dpng','Binary_Sparse_Image_After_Hole_Filling.png');
    end
else
    E6 = im2bw(E5);
end

objCnt = 1; Mask = [];
for i=1:length(MaskIdxList),
    Temp_Mask = zeros(Row, Col);
    % Copy only the masked pixels in
    Temp_Mask(MaskIdxList(i).PixelIdxList) = 1;
    if sum(Temp_Mask(:)==1) >= connectivity_area && sum(Temp_Mask(:)==1) < 0.8*Row*Col,
        Mask{objCnt} = Temp_Mask;
        objCnt = objCnt + 1;
    end
end

if isempty(Mask),
    disp('No object found!');
    return;
end

if IsDebugging
	% for ii = 1:length(bxFinal)
	for ii = 1:length(MaskIdxList) % Added on 2015/11/07
		figure;imshow(Mask{ii},[]);title(strcat('SubFrame Mask for object', num2str(ii))); size(Mask{ii})
		print(gcf,'-dpng',strcat('SubFrame_Mask_',num2str(ii),'.png'));
	end
end