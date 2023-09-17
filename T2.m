
function [T2_Map ] = T2(varargin)

%get all images
IM1 = dicomread('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0138');
IM2 = dicomread('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0139');
IM3 = dicomread('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0140');
IM4 = dicomread('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0141');
IM5 = dicomread('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0142');

header = dicominfo('/Users/neuroimaging/Desktop/Imaging 1 Practicum/Sammet_Lab/Lab_4_Scans/Volunteer 5/DICOM/IM_0138');
nx     = header.Height;
ny     = header.Width;

%the most inelegant code, but it'll get it done hopefully

%pixel rebinning for time saving....
%nx = 224;
%ny = 224;
Total_Slices = 5;
Image_stack=zeros(Total_Slices, nx, ny);
%stack images
Image_stack(1,:,:) = imresize(IM1,[nx ny],'bicubic');
Image_stack(2,:,:) = imresize(IM2,[nx ny],'bicubic');
Image_stack(3,:,:) = imresize(IM3,[nx ny],'bicubic');
Image_stack(4,:,:) = imresize(IM4,[nx ny],'bicubic');
Image_stack(5,:,:) = imresize(IM5,[nx ny],'bicubic');

T2_Map = zeros(nx,ny);
%now fitting to T2 decay... 
%for pixels in rectangular FOV chosen by eye
for i = 99:345
    for j = 84:367
        if Image_stack(1, i, j)>0
            T2_Signal = double(Image_stack(1:Total_Slices, i,j)); %get signal from particular voxel for all images along z axis
            [T2_val, f0, TEs, Signal] = T2fit(T2_Signal);
            %plot(f0, TEs, Signal./Signal(1))
            %title(strcat(string(i),string(j)));
            T2_Map(i,j) = T2_val;
        else
            T2_Map(i,j)=0;
        end
    end
    counter = strcat(string(i),string(j)) %count/keep track of which voxel we're at.
end

figure,imshow(T2_Map,[])


function [T2_val, f0, TEs, Signal] = T2fit(Signal)
    TEs = double([20, 40, 60, 80, 100]); %time in ms (as each TE is 20 ms??)
    startpoints = [1, 90]; %starting unnormalized value of 600 or so, and a T2 of bone or cartilage of 90?
    T2_equation = 'a*exp(-x/b)'; %normalized.... x = TE, b = T2, a = whatever constant
    [f0, G0] = fit(TEs', Signal./Signal(1), T2_equation, 'Start', startpoints); %fit normalized!
    T2_val = f0.b;
end

end