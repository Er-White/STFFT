clear; clc;
maindir_SM_downscale = dir('G:\PHD\18-19neimeng\SM_downscale\*.tif');
maindir_SMAP = dir('G:\PHD\18-19neimeng\SMAP_CLIP\*.tif');
maindir_output = "G:\PHD\18-19neimeng\FUSION1\";
n = 10 %低频数量
Fs = 1;            % Sampling frequency                    
T = 1/Fs;          % Sampling period       
L = 183;             % Length of signal
t = (0:L-1)*T;        % Time vector
for i = 184:366
    subdir_SMAP = fullfile(maindir_SMAP(i).folder, maindir_SMAP(i).name);
    subdir_SMAP = convertCharsToStrings(subdir_SMAP);
    [SMAP,R1] = readgeoraster(subdir_SMAP,"OutputType","double");
    data_train_temp = reshape(SMAP,1,[]);
    data_train(i-183,:) = data_train_temp;
end
[m1,n1] = size(SMAP);
for i = 1:n1*m1
    sig_fft = fft(data_train(:,i));
    sig_fft=sig_fft/L;   % 将频域序列 X 除以序列的长度 N
    P1 = sig_fft(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    A(:,i)=abs(P1);           % 计算频域序列 Y 的幅值
    Pha(:,i)=angle(P1);       % 计算频域序列 Y 的相角 (弧度制)
    f = Fs*(0:(L/2))/L;
    f =f';
end
A = A'; Pha = Pha'; 
A = reshape(A,m1,n1,[]); Pha = reshape(Pha,m1,n1,[]);

COL = [];
for i = 15:27
    subdir_SM_downscale = fullfile(maindir_SM_downscale(i).folder, maindir_SM_downscale(i).name);
    subdir_SM_downscale = convertCharsToStrings(subdir_SM_downscale);
    [SM_downscale,R2] = readgeoraster(subdir_SM_downscale,"OutputType","double");
    COL = cat(3,COL,SM_downscale);
end
[m2,n2] = size(SM_downscale);
A = imresize(A, [m2,n2], "bicubic"); 
F = f;
Pha = mean(Pha,[1 2]);
Pha = reshape(Pha, [] ,1);
% Pha = imresize(Pha, [m2,n2], "nearest");
A(:,:,1:n) = 1;
% time = [3 15 27 39 51 63 75 87 99 111 123 135 159 183];
time = [10 22 34 46 58 70 82 94 106 118 130 154 178];
for i = 1:n1*m1
    if rem(i,m1) == 0
        rowRange = (m1 - 1) * 450 + 1 : m1 * 450;
    else
        rowRange = (rem(i,m1) - 1) * 450 + 1 : rem(i,m1) * 450;
    end
    colRange = (ceil(i/m1) - 1) * 450 + 1 : ceil(i/m1) * 450;
    COL_temp = COL(rowRange, colRange,:);
    A_temp = A(rowRange, colRange, n+1:end);
    A_temp = reshape(A_temp,[],length(f)-n);
    result_c = A_temp*cos(2*pi*F(n+1:end)*t+Pha(n+1:end));
    result_c = result_c(:,time);
    COL_temp = reshape(COL_temp,[],length(time))-result_c;
    COL_temp = COL_temp';
    result_half = cos(2*pi*F(1:n)*t(time)+Pha(1:n))';
    C = result_half\COL_temp;
    C = [C',A_temp];
    result = C*cos(2*pi*F*t+Pha);
    result = reshape(result,450,450,183);
    SMC(rowRange, colRange,:) = result;
end

for i = 184:366
    output_name = maindir_output + maindir_SMAP(i).name;
    geotiffwrite(output_name, SMC(:,:,i-183), R2);
end
SMC_estimated_t = reshape(SMC(:,:,time(1)),[],1);
ccc = reshape(COL(:,:,1),[],1);
scatter(SMC_estimated_t,ccc)
latitude = 40.869177; % Replace with your desired latitude
longitude = 107.129053; % Replace with your desired longitude
[xIntrinsic,yIntrinsic] = geographicToIntrinsic(R2,latitude,longitude);
index = [xIntrinsic,yIntrinsic];
index = round(index);
data = SMC(index(2),index(1),1:183);
data = reshape(data,1,[]);
data = data';
plot(data)
