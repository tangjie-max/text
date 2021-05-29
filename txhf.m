% ----------------2����һ���Ҷ�ͼ������˶�ģ�������Ӹ�˹������������ά���˲����и�ԭ-------------
clc;         %��տ���̨
clear;       %��չ�����
close all;   %�ر��Ѵ򿪵�figureͼ�񴰿�
color_pic=imread('zx.jpg');  %��ȡ��ɫͼ��
gray_pic=rgb2gray(color_pic);    %����ɫͼת���ɻҶ�ͼ
double_gray_pic=im2double(gray_pic);   %��uint8ת��im2double�ͱ��ں��ڼ���
[width,height]=size(double_gray_pic);

%-------------------------����˶�ģ��----------------------
H_motion = fspecial('motion', 18, 90);%�˶�����Ϊ18����ʱ���˶��Ƕ�Ϊ90��
motion_blur = imfilter(double_gray_pic, H_motion, 'conv', 'circular');%����˲�
noise_mean=0;  %��Ӿ�ֵΪ0
noise_var=0.001; %����Ϊ0.001�ĸ�˹����
motion_blur_noise=imnoise(motion_blur,'gaussian',noise_mean,noise_var);%��Ӿ�ֵΪ0������Ϊ0.001�ĸ�˹����
figure('name','�˶�ģ������');
subplot(1,2,1);imshow(motion_blur,[]);title('�˶�ģ��');
subplot(1,2,2);imshow(motion_blur_noise,[]);title('�˶�ģ���������');

%----------------------------ά���˲�matlab�Դ�����deconvwnr-------------------------
restore_ignore_noise = deconvwnr(motion_blur_noise, H_motion, 0);  %nsr=0,��������
signal_var=var(double_gray_pic(:));
estimate_nsr=noise_var/signal_var;   %���űȹ�ֵ
restore_with_noise=deconvwnr(motion_blur_noise,H_motion,estimate_nsr);  %�źŵĹ�����ʹ��ͼ��ķ�����ƹ���
figure('name','������ά���˲�');
subplot(1,2,1);imshow(im2uint8(restore_ignore_noise),[]);title('��������ֱ��ά���˲�(nsr=0)���൱�����˲�');
subplot(1,2,2);imshow(im2uint8(restore_with_noise),[]);title('��������ά���˲�');

%---------------------------��ʽ��----------------------------
fourier_H_motion=fft2(H_motion,width,height);  %H(u,v)
pow_H_motion=abs(fourier_H_motion).^2;   %|H(u,v)|^2
noise=motion_blur_noise-motion_blur;    %��ȡ��������
fourier_noise=fft2(noise);    % N(u,v)  ��������Ҷ�任
fourier_double_gray_pic=fft2(double_gray_pic);  %F(u,v)Ϊδ�����˻���ͼƬ
nsr=abs(fourier_noise).^2./abs(fourier_double_gray_pic).^2;   %���ű�=|N(u,v)|^2/|F(u,v)|^2
H_w=1./fourier_H_motion.*pow_H_motion./(pow_H_motion+nsr); %H_w(u,v)=1/H(u,v)*|H(u,v)|^2/[|H(u,v)|^2+NSR] 
fourier_motion_blur_noise=fft2(motion_blur_noise);  %G(u,v)
restore_with_noise=ifft2(fourier_motion_blur_noise.*H_w);  %���Ƶ��=G(u,v)H_w(u,v)��ʱ��ΪƵ����Ҷ��任
figure('name','��ʽ��άҲ���˲�');
subplot(1,2,1);imshow(motion_blur_noise,[]);title('�˶�ģ���������')
subplot(1,2,2);imshow(restore_with_noise,[]);title('άҲ���˲�')

