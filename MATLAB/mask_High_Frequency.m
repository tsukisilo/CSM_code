%低通滤波mask
function filted_img_Frequency=mask_High_Frequency(img_Frequency)
[rows,cols]=size(img_Frequency);
radius=min(rows,cols)*0.25;
center_x=(cols+1)/2;%列坐标
center_y=(rows+1)/2;%行坐标
for i=1:rows
    for j=1:cols
        if(sqrt((i-center_y)^2+(j-center_x)^2)>radius)
            img_Frequency(i,j)=0;
        end
    end
end
filted_img_Frequency=img_Frequency;