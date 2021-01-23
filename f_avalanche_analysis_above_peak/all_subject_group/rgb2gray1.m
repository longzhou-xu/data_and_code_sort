function rgb2gray1(FigPath)
%
%
MyOriginalFig = imread([FigPath,'.tif']);
MyGrayFig = rgb2gray(MyOriginalFig);
imwrite(MyGrayFig , [FigPath,'newgray.tif'], 'tif');  
  
%显示原来的RGB图像  
figure(1);  
imshow(MyOriginalFig);  
  
%显示经过系统函数运算过的灰度图像  
figure(2);  
imshow(MyGrayFig);
end

