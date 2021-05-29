% ----------------2、对一幅灰度图像进行运动模糊并叠加高斯噪声，并采用维纳滤波进行复原-------------
clc;         %清空控制台
clear;       %清空工作区
close all;   %关闭已打开的figure图像窗口
color_pic=imread('zx.jpg');  %读取彩色图像
gray_pic=rgb2gray(color_pic);    %将彩色图转换成灰度图
double_gray_pic=im2double(gray_pic);   %将uint8转成im2double型便于后期计算
[width,height]=size(double_gray_pic);

%-------------------------添加运动模糊----------------------
H_motion = fspecial('motion', 18, 90);%运动长度为18，逆时针运动角度为90°
motion_blur = imfilter(double_gray_pic, H_motion, 'conv', 'circular');%卷积滤波
noise_mean=0;  %添加均值为0
noise_var=0.001; %方差为0.001的高斯噪声
motion_blur_noise=imnoise(motion_blur,'gaussian',noise_mean,noise_var);%添加均值为0，方差为0.001的高斯噪声
figure('name','运动模糊加噪');
subplot(1,2,1);imshow(motion_blur,[]);title('运动模糊');
subplot(1,2,2);imshow(motion_blur_noise,[]);title('运动模糊添加噪声');

%----------------------------维纳滤波matlab自带函数deconvwnr-------------------------
restore_ignore_noise = deconvwnr(motion_blur_noise, H_motion, 0);  %nsr=0,忽视噪声
signal_var=var(double_gray_pic(:));
estimate_nsr=noise_var/signal_var;   %噪信比估值
restore_with_noise=deconvwnr(motion_blur_noise,H_motion,estimate_nsr);  %信号的功率谱使用图像的方差近似估计
figure('name','函数法维纳滤波');
subplot(1,2,1);imshow(im2uint8(restore_ignore_noise),[]);title('忽视噪声直接维纳滤波(nsr=0)，相当于逆滤波');
subplot(1,2,2);imshow(im2uint8(restore_with_noise),[]);title('考虑噪声维纳滤波');

%---------------------------公式法----------------------------
fourier_H_motion=fft2(H_motion,width,height);  %H(u,v)
pow_H_motion=abs(fourier_H_motion).^2;   %|H(u,v)|^2
noise=motion_blur_noise-motion_blur;    %提取噪声分量
fourier_noise=fft2(noise);    % N(u,v)  噪声傅里叶变换
fourier_double_gray_pic=fft2(double_gray_pic);  %F(u,v)为未经过退化的图片
nsr=abs(fourier_noise).^2./abs(fourier_double_gray_pic).^2;   %噪信比=|N(u,v)|^2/|F(u,v)|^2
H_w=1./fourier_H_motion.*pow_H_motion./(pow_H_motion+nsr); %H_w(u,v)=1/H(u,v)*|H(u,v)|^2/[|H(u,v)|^2+NSR] 
fourier_motion_blur_noise=fft2(motion_blur_noise);  %G(u,v)
restore_with_noise=ifft2(fourier_motion_blur_noise.*H_w);  %输出频域=G(u,v)H_w(u,v)，时域为频域傅里叶逆变换
figure('name','公式法维也纳滤波');
subplot(1,2,1);imshow(motion_blur_noise,[]);title('运动模糊添加噪声')
subplot(1,2,2);imshow(restore_with_noise,[]);title('维也纳滤波')

