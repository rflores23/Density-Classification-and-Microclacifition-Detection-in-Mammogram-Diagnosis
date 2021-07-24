clear; clc; close all

fprintf("IMAGE ACQUISITION\n")
cd '/Users/rickflores/Documents/BI 260/mdb200data'
pgmFiles = dir('/Users/rickflores/Documents/BI 260/mdb200data/*.pgm');



% I = imread('mdb188.pgm');

imgSize = 1024;
num = 200;
imgList = zeros(imgSize,imgSize,num);
for i = 1:num
    fprintf("reading image " + i + "\n")
    imgList(:,:,i) = imread(pgmFiles(i).name);
end
imgList = double(imgList) ./ max(imgList,[],'all');

I = [];
Tbl = zeros(num,3);
for i = 1:num
    I = imgList(:,:,i);
    I = I(200:900,350:900);
S = strel('square', 3);
phi = imclose(I, S);
gamma = imopen(phi, S);


Min = min(gamma, I);

T = I - Min;
T_double = im2double(T);
calc = im2bw(T_double, (4/255));
% figure()
% subplot(132)
% imshow(calc)

%erode away small spots
J = imerode(calc, strel('disk', 1));
fin = imreconstruct(J, calc);
% subplot(133)
% imshow(fin)
% 
% 
% subplot(131)
% imshow(I)

%find number of calcifications
CC = bwconncomp(fin);
num_calcifications = CC.NumObjects;

%find average microcalcification size
objects = CC.PixelIdxList;
size_objects = [];
for i = 1:length(objects)
    size_objects = [size_objects, size(objects{i})];
end
average_calc_size = mean(size_objects);

%find avergae distance between calcification
object_locations = [];
%first center location of each object
for i = 1:length(objects)
        object_locations = [object_locations, mean(objects{i})];
end

clustering = norm(object_locations)./10e4;

Tbl(i,:) = [num_calcifications, average_calc_size, clustering];

end