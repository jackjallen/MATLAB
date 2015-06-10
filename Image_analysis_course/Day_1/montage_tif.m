function montage_tif(filename, starting_slice)

%filename = 'mri-stack.tif'
image_info = imfinfo(filename)

image_width = image_info.Width;
image_height = image_info.Height;

image_number_of_stacks = size(filename);
preallocated_matrix = zeros(image_height, image_width, image_number_of_stacks(1));

for i = 1:image_number_of_stacks
   
    preallocated_matrix(:,:,i) = imread(filename);
    
end

imagesc(preallocated_matrix(:,:,1))

montage(preallocated_matrix./image_number_of_stacks(1), 1, 'indices', starting_slice:image_number_of_stacks) 

end