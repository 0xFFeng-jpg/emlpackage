function Pitch = getPitch(Img)
if(1)
figure;imshow(Img,[]);
rect = getrect;
Pitch = Img(rect(2):rect(2)+rect(4),rect(1):rect(1)+rect(3),:);
else
    Pitch = Img(6:282,54:286,:);
end


