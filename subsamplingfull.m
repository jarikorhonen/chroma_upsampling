function psnr = subsamplingfull(imfile, Q)

% Input: imfile: filename for the input image
%        Q: quantization factor for simulated JPEG quantization
%           valid values [1,100], otherwise no JPEG simulated
%
% Output: vector of psnr values
%         [proposed box bilinear bicubic laczos3]  
%
%

% Read image from file and normalize to [0,1]
ref_image = cast(imread(imfile),'double')./255;

% Downsample in RGB domain to avoid chroma loss impact in original image
%ref_image = ref_image(1:end/4,1:end/4,:);
%ref_image = imresize(ref_image, 0.5);
%ref_image = (round(ref_image .* 255.0))./255.0;

% Here simulate JPEG compression (additional)
% ref_image(:,:,1) = luma_quant(ref_image(:,:,1),50);
% ref_image(:,:,2) = luma_quant(ref_image(:,:,2),50);
% ref_image(:,:,3) = luma_quant(ref_image(:,:,3),50);
% ref_image = (round(ref_image .* 255.0))./255.0;


% Make sure that the resolution is divisible by two
dim = size(ref_image(:,:,1));
if mod(dim(1),2)==1
    ref_image = ref_image(1:end-1,:,:);
end
if mod(dim(2),2)==1
    ref_image = ref_image(:,1:end-1,:);
end

% Covert from RGB to Y'CbCr 
ycbcr_image = rgb2ycbcr(ref_image);

y = ycbcr_image(:,:,1);
cb = downsample(ycbcr_image(:,:,2));
cr = downsample(ycbcr_image(:,:,3));

% Here simulate JPEG compression (optional)
if Q>0 && Q<101
    y = luma_quant(y,Q);
    cb = chroma_quant(cb,Q);
    cr = chroma_quant(cr,Q);
else
    disp('Invalid Q, no JPEG compression simulated');
end

psnr = [];
ssim = [];

% Perform proposed upsampling and compute psnr
ycbcr_im_ups = upsample(y,cb,cr);
rgb_image = ycbcr2rgb(ycbcr_im_ups);
rgb_image = (round(rgb_image .* 255.0))./255.0;
psnr(1) = 10*log10(1/mean(mean(mean((rgb_image-ref_image).^2))));
psnr(2) = SSIMfig(rgb_image*255, ref_image*255);

%figure, imshow(ref_image);
%figure, imshow(rgb_image);


% psnr(2): box
ycbcr_im_ups(:,:,2) = imresize(cb,2.0,'box'); 
ycbcr_im_ups(:,:,3) = imresize(cr,2.0,'box'); 
rgb_image = ycbcr2rgb(ycbcr_im_ups);
rgb_image = (round(rgb_image .* 255.0))./255.0;
psnr(3) = 10*log10(1/mean(mean(mean((rgb_image-ref_image).^2))));
psnr(4) = SSIMfig(rgb_image*255, ref_image*255);

% psnr(3): bilinear
ycbcr_im_ups(:,:,2) = imresize(cb,2.0,'bilinear'); 
ycbcr_im_ups(:,:,3) = imresize(cr,2.0,'bilinear'); 
rgb_image = ycbcr2rgb(ycbcr_im_ups);
rgb_image = (round(rgb_image .* 255.0))./255.0;
psnr(5) = 10*log10(1/mean(mean(mean((rgb_image-ref_image).^2))));
psnr(6) = SSIMfig(rgb_image*255, ref_image*255);

% psnr(4): bicubic
ycbcr_im_ups(:,:,2) = imresize(cb,2.0,'bicubic'); 
ycbcr_im_ups(:,:,3) = imresize(cr,2.0,'bicubic'); 
rgb_image = ycbcr2rgb(ycbcr_im_ups);
rgb_image = (round(rgb_image .* 255.0))./255.0;
psnr(7) = 10*log10(1/mean(mean(mean((rgb_image-ref_image).^2))));
psnr(8) = SSIMfig(rgb_image*255, ref_image*255);

% psnr(5): Lanczos-3
ycbcr_im_ups(:,:,2) = imresize(cb,2.0,'lanczos3'); 
ycbcr_im_ups(:,:,3) = imresize(cr,2.0,'lanczos3'); 
rgb_image = ycbcr2rgb(ycbcr_im_ups);
rgb_image = (round(rgb_image .* 255.0))./255.0;
psnr(9) = 10*log10(1/mean(mean(mean((rgb_image-ref_image).^2))));
psnr(10) = SSIMfig(rgb_image*255, ref_image*255);

return;

function img = luma_quant(in_y, Q)

% Simulate JPEG compression for luma with quantization factor Q
QTable_Y = [16 11 10 16  24  40  51  51
            12 12 14 19  26  58  60  55
            14 13 16 24  40  57  69  56
            14 17 22 29  51  87  80  62
            18 22 37 56  68 109 103  77
            24 35 55 64  81 104 113  92
            49 64 78 87 103 121 120 101
            72 92 95 98 112 100 103  99];

if Q < 50
    Q = floor(5000/Q);
else
    Q = 200-2*Q;
end
QTable_Y = round(cast(QTable_Y,'double').*Q/100.0+0.5);
        
dim = size(in_y);

% This is to avoid problems in case dimensions not divisable by 8
img = in_y;

for i=1:8:dim(1)-7
    for j=1:8:dim(2)-7
        QBlock = dct2(in_y(i:i+7,j:j+7).*255);   % forward transform
        QBlock = floor(QBlock./(QTable_Y)+0.5);  % quantization 
        %QBlock = round(16*QBlock/Q./QTable_Y);
        iQBlock = QBlock.*QTable_Y;              % inverse quantization 
        %iQBlock = (QBlock+0.5*sign(QBlock)).*QTable_Y*Q/16;
        iQBlock = idct2(iQBlock);                % inverse transform
        img(i:i+7,j:j+7)=iQBlock./255;
    end
end

return;

function img = chroma_quant(in_c, Q)

% Simulate JPEG compression for luma with quantization factor Q
QTable_C = [17 18 24 47 66 99 99 99
            18 21 26 66 99 99 99 99
            24 26 56 99 99 99 99 99
            47 66 99 99 99 99 99 99
            99 99 99 99 99 99 99 99
            99 99 99 99 99 99 99 99
            99 99 99 99 99 99 99 99
            99 99 99 99 99 99 99 99];

if Q < 50
    Q = floor(5000/Q);
else
    Q = 200-2*Q;
end
QTable_C = round(cast(QTable_C,'double').*Q/100.0+0.5);        
        
dim = size(in_c);

% This is to avoid problems in case dimensions not divisable by 8
img = in_c;

for i=1:8:dim(1)-7
    for j=1:8:dim(2)-7
        QBlock = dct2(in_c(i:i+7,j:j+7).*255);   % forward transform
        QBlock = floor(QBlock./(QTable_C)+0.5);  % quantization 
        iQBlock = QBlock.*QTable_C;              % inverse quantization 
        iQBlock = idct2(iQBlock);                % inverse transform
        img(i:i+7,j:j+7)=iQBlock./255;
    end
end

return;

function img = downsample(in_c, in_y)

% Simple chroma downscaling by averaging from values in 2x2 blocks
dim = size(in_c);
img = zeros(dim(1)/2, dim(2)/2);

for i=1:2:dim(1)-1
    for j=1:2:dim(2)-1
        img((i-1)/2+1,(j-1)/2+1) = ...
            (in_c(i,j)+in_c(i,j+1)+in_c(i+1,j)+in_c(i+1,j+1))/4.0;
    end
end

return;

function img = upsample(in_y, in_cb, in_cr)

% This function implements the proposed upsampling scheme

dim = size(in_y);

% Luma values will be returned unchanged
img(:,:,1) = in_y;

% Bilinear interpolation needed to generate pixels at edges
img(:,:,2) = imresize(in_cb,2,'bilinear');
img(:,:,3) = imresize(in_cr,2,'bilinear');

totdiff = 0;
n = 0;

% Go through all the (luma) pixels, except edge pixels
for i=2:dim(1)-1
    for j=2:dim(2)-1
        
        weight=[0 0 0 0]; % Weights for four blocks
        msd=[0 0 0 0];    % MSDs for four blocks
        c_pos=zeros(4,2);  % Chroma sample positions for the four blocks
        
        if mod(i,2)==0 && mod(j,2)==0
            % Lower right pixel in 2x2 block: nearest chromas to right,
            % to down, and to downright
            msd(1) = mean(mean((in_y(i,j)-in_y(i-1:i,   j-1:j)).^2));
            msd(2) = mean(mean((in_y(i,j)-in_y(i+1:i+2, j-1:j)).^2));
            msd(3) = mean(mean((in_y(i,j)-in_y(i-1:i,   j+1:j+2)).^2));
            msd(4) = mean(mean((in_y(i,j)-in_y(i+1:i+2, j+1:j+2)).^2));
            c_pos(1,:) = [i/2   j/2];
            c_pos(2,:) = [i/2+1 j/2];
            c_pos(3,:) = [i/2   j/2+1];
            c_pos(4,:) = [i/2+1 j/2+1];

        elseif mod(i,2)==1 && mod(j,2)==0
            
            % Upper right pixel in 2x2 block: nearest chromas to right,
            % to up, and to upright
            msd(1) = mean(mean((in_y(i,j)-in_y(i:i+1,   j-1:j)).^2));
            msd(2) = mean(mean((in_y(i,j)-in_y(i-2:i-1, j-1:j)).^2));
            msd(3) = mean(mean((in_y(i,j)-in_y(i:i+1,   j+1:j+2)).^2));
            msd(4) = mean(mean((in_y(i,j)-in_y(i-2:i-1, j+1:j+2)).^2));
            c_pos(1,:) = [ceil(i/2)  j/2];
            c_pos(2,:) = [floor(i/2) j/2];
            c_pos(3,:) = [ceil(i/2)  j/2+1];
            c_pos(4,:) = [floor(i/2) j/2+1];

        elseif mod(i,2)==0 && mod(j,2)==1
            % Lower left pixel in 2x2 block: nearest chromas to left,
            % to down, and to downleft
            msd(1) = mean(mean((in_y(i,j)-in_y(i-1:i,   j:j+1)).^2));
            msd(2) = mean(mean((in_y(i,j)-in_y(i+1:i+2, j:j+1)).^2));
            msd(3) = mean(mean((in_y(i,j)-in_y(i-1:i,   j-2:j-1)).^2));
            msd(4) = mean(mean((in_y(i,j)-in_y(i+1:i+2, j-2:j-1)).^2));
            c_pos(1,:) = [i/2   ceil(j/2)];
            c_pos(2,:) = [i/2+1 ceil(j/2)];
            c_pos(3,:) = [i/2   floor(j/2)];
            c_pos(4,:) = [i/2+1 floor(j/2)];

        else
            % Upper left pixel in 2x2 block: nearest chromas to left,
            % to up, and to upleft            
            msd(1) = mean(mean((in_y(i,j)-in_y(i:i+1,   j:j+1)).^2));
            msd(2) = mean(mean((in_y(i,j)-in_y(i-2:i-1, j:j+1)).^2));
            msd(3) = mean(mean((in_y(i,j)-in_y(i:i+1,   j-2:j-1)).^2));
            msd(4) = mean(mean((in_y(i,j)-in_y(i-2:i-1, j-2:j-1)).^2));
            c_pos(1,:) = [ceil(i/2)  ceil(j/2)];
            c_pos(2,:) = [floor(i/2) ceil(j/2)];
            c_pos(3,:) = [ceil(i/2)  floor(j/2)];
            c_pos(4,:) = [floor(i/2) floor(j/2)];
           
        end
        
        % Compute exponent alfa
        alfa = 1;
        theta = 0.15;
        diff2 = 0;
        if max(msd)-min(msd)<theta
            alfa=(max(msd)-min(msd))/theta;
        end

        % Set initial weights from MSD values
        [d ind] = sort(msd,'descend');
        weight(ind(1))=msd(ind(4));
        weight(ind(2))=msd(ind(3));
        weight(ind(3))=msd(ind(2));
        weight(ind(4))=msd(ind(1));
        
        % Combine distance weights and MSD weights
        weight = weight.^alfa;
        weight = weight.*[0.5625 0.1875 0.1875 0.0625].^(1-alfa);

        % Normalize
        weight = weight./sum(weight);

        % Create upsampled Cb pixel
        pix(1) = in_cb(c_pos(1,1),c_pos(1,2)); 
        pix(2) = in_cb(c_pos(2,1),c_pos(2,2));
        pix(3) = in_cb(c_pos(3,1),c_pos(3,2)); 
        pix(4) = in_cb(c_pos(4,1),c_pos(4,2));
        img(i,j,2)=sum(weight.*pix);
        
        % Create upsampled Cr pixel
        pix(1) = in_cr(c_pos(1,1),c_pos(1,2)); 
        pix(2) = in_cr(c_pos(2,1),c_pos(2,2));
        pix(3) = in_cr(c_pos(3,1),c_pos(3,2)); 
        pix(4) = in_cr(c_pos(4,1),c_pos(4,2));
        img(i,j,3)=sum(weight.*pix);
    end
end

return;