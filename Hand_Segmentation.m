% Yusuf Gökberk Keptiğ 2079366

clear all;
clc;
close all;

% % Reading all files from dataset
 db= imageDatastore('\Dataset\','IncludeSubfolders',true,'LabelSource','foldernames','FileExtensions', {'.png'});
 images=readall(db);


for n=1:20
[s1,s2,s3]=size(images{n});
% figure; imshow(I);

R = double(images{n}(:,:,1));
G = double(images{n}(:,:,2));
B = double(images{n}(:,:,3));

r=  R./(R+G+B);
b=  B./(R+G+B);
g = G./(R+G+B);

% taking R-plane, B-plane, G-plane values as features
X = [r(:) g(:) b(:)];

r2=uint8(r.*255);
b2=uint8(b.*255);
g2=uint8(g.*255);
y{n}(:,:,1)=r2;
y{n}(:,:,2)=g2;
y{n}(:,:,3)=b2;

yR = double(y{n}(:,:,1));
yG = double(y{n}(:,:,2));
yB = double(y{n}(:,:,3));

%Get histValues for each channel
[yRed, valueRed] = imhist(r2);
[yGreen, valueGreen] = imhist(g2);
[yBlue, valueBlue] = imhist(b2);


%Plot them together in one plot
% plot(valueRed, yRed, 'Red', valueGreen, yGreen, 'Green', valueBlue, yBlue, 'Blue');



%% Finding Peaks 
% we take the peaks then we sort them in descending order after that we
% find index that corresponds each 3 largest value using findValues
% function

pksB = findpeaks(yBlue);
sPksB = sort(pksB,'descend');
blueP1 = sPksB(1);
blueP2 = sPksB(2);
blueP3 = sPksB(3);
[bluePeak1,bluePeak2,bluePeak3]=findValues(yBlue,blueP1,blueP2,blueP3);


pksR = findpeaks(yRed);
sPksR = sort(pksR,'descend');
redP1 = sPksR(1);
redP2 = sPksR(2);
redP3 = sPksR(3);
[redPeak1,redPeak2,redPeak3]=findValues(yRed,redP1,redP2,redP3);


pksG = findpeaks(yGreen);
sPksG = sort(pksG,'descend');
greenP1 = sPksG(1);
greenP2 = sPksG(2);
% sometimes green doesn't have 3 peaks so this is an alternative solution
% for that
if(size(sPksG)>2)
   greenP3 = sPksG(3); 
else
     greenP3 =greenP2;
end

[greenPeak1,greenPeak2,greenPeak3]=findValues(yGreen,greenP1,greenP2,greenP3);


centroid1 = [redPeak1,greenPeak3,bluePeak3];
centroid2 = [redPeak2,greenPeak2,bluePeak3];
centroid3 = [redPeak3,greenPeak3,bluePeak3];

%% image writing 

clusterGroup = Kmeans(centroid1,centroid2,centroid3,yR,yG,yB,s1,s2);

outimg = Masking(clusterGroup,s1,s2);

outputimg = uint8(outimg);
filename = sprintf('mask_image%d.png',n);
FILENAME = ['\Dataset\SegmentationResults\', filename];
imwrite(outputimg, FILENAME);


end




%%  findValues Function 
% this function finds 3 max peak values' indexis
% traversing arrray and finding index that are corresponding to their
% values 

function [first,second,third] = findValues(arr,f,s,th)
    first=0;
    second=0;
    third=0;
    for i=1:255
       if(arr(i)==f)
           first = i-1;
       end
       if(arr(i)==s)
           second = i-1;
       end
       if(arr(i)==th)
           third = i-1;
       end
        
    end
    
end

%% Masking Funtion
function [outimg] = Masking(clusterG,s1,s2)

outimg=zeros(s1,s2);
    for i=1:s1
        for j=1:s2
            if clusterG(i,j)== 1
                outimg(i,j)= 62;
            elseif clusterG(i,j)== 2
                outimg(i,j)= 124;
            elseif clusterG(i,j)== 3
                outimg(i,j)= 0;
            end
        end
    end
figure;imshow(uint8(outimg)); % segmented image

end


%% Kmeans Function
function [clusterGroup]=Kmeans(centroid1,centroid2,centroid3,yR,yG,yB,s1,s2)

temp1= centroid1;
temp2= centroid2;
temp3= centroid3;

c=0;


clusterGroup=zeros(s1,s2);
group1=0;
group2=0;
group3=0;
cluster1{1}=centroid1; 
cluster2{1}=centroid2; 
cluster3{1}=centroid3; 

while(c==0)
    c1=temp1;
    c2=temp2;
    c3=temp3;
    for i=1:s1
        for j=1:s2
            pVector = [yR(i,j) ,yG(i,j) ,yB(i,j)];
            
            %Euclidean distance between rgb values 
            if(norm(temp1-pVector)>norm(temp2-pVector) && norm(temp1-pVector)>norm(temp3-pVector) )
                clusterGroup(i,j) = 1;
                group1 = group1+1;
                cluster1{group1+1}= pVector;
                
            elseif(norm(temp2-pVector)>norm(temp1-pVector) && norm(temp2-pVector)>norm(temp3-pVector) )
                clusterGroup(i,j) = 2;
                group2 = group2+1;
                cluster2{group2+1}= pVector;
                 
            elseif(norm(temp3-pVector)>norm(temp1-pVector) && norm(temp3-pVector)>norm(temp1-pVector) )
                clusterGroup(i,j) = 3;
                group3 = group3+1;
                cluster3{group3+1}= pVector;
            end
        end
    end
%     in order to use mean Function we need to change cell2mat then find
%     mean
    t1C = cell2mat(cluster1);
    temp1 = mean(t1C);
    t2C = cell2mat(cluster2);
    temp2 = mean(t2C);
    t3C = cell2mat(cluster3);
    temp3 = mean(t3C);

    
%     if the difference is less than threshold value continues otherwise
%     stops the iteration
    if(norm(temp1-c1)<= 0.5 && norm(temp2-c2) <= 0.5 && norm(temp3-c3)<= 0.5)
       c1={};
       c2={};
       c3={};
       
    else
       c = 1;
    end
end

end
