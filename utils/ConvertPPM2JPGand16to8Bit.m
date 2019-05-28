function ConvertPPM2JPGand16to8Bit
% Original images in the database had to be converted to 8bit grayscale
% images. Original format was PPM 16bit.
% 
% Used only once at the beginning of the project.
%
%
directory = 'E:\\Eye movement data for Whitening study\\Natural Images database\\';

listing = dir(directory);

usefulFiles = zeros(size(listing));
for i=1:length(listing)
    if ~(listing(i).isdir)
        [~,~,ext] = fileparts(listing(i).name);
        if ~isempty(strfind(ext,'ppm'))
            usefulFiles(i) = 1;
            
            % read the image
            fullfilename = [directory listing(i).name];
            im = imread(fullfilename);
            
            % convert it to 8bit
            im2 = uint8(im/256);
            clear im;
            
            % gamma correction
            im3 = imadjust(im2,[],[],0.5);
            clear im2;
            
            % rgb to grayscale conversion
            im4 = rgb2gray(im3);
            clear im3;
            
            % save as jpg
            newFileName = [fullfilename(1:end-4) '.jpg'];
            imwrite(im4,newFileName);
            
        end
    end
end