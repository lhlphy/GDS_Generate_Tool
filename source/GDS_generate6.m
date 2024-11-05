%Code Author:Li Haolin
%Last revised date:2024/11/5
%clear;
%单区超表面，几何相位和传播相位共同作用，左旋和右旋偏振光全息图像不同，全部显示绿光
% set the work path to parentFolder
currentFile = mfilename('fullpath');  % 获取当前脚本或文件的完整路径
currentFolder = fileparts(currentFile); % 获取文件所在文件夹的路径
parentFolder = fileparts(currentFolder);  % 获取上一层文件夹的路径
cd([parentFolder,'\Generated']); % 设置工作目录为上一层文件夹

[log,cad,cellname,topcell,layerMap] = InitializeCell();
refs = [];
refs = GetRefsFloorplan(refs);


Px = 1; %大单元格x周期
Py = 1; %大单元格y周期
a = 0.2*[0,2.35,1.85,1.40]; %矩形块长轴
b = 0.06*[0,2.35,1.85,1.40]; %矩形块短轴
Rotation = Comb_hologram*180/(2*pi)/2; %注意几何相位为转动角Rot的两倍
lam_eff = 0.2 ; %SPP波长 %请计算
Phi_grad = pi/2; %相位梯度 rad/um
lam_SPP =  [100,LS(620*1e-9,1),LS(530*1e-9,1),LS(470*1e-9,1)]; %SPP波长 %请计算
% imgmax = max(max(Img_color)); %色彩匹配编号
% Img_color = Img_color +(4-imgmax);  %将0，1，2，3转为1，2，3，4
Img_C = 3*ones(size(Img_color));

cell1 = gds_structure('cell1');
%开始逐个单元绘制版图
[NX,NY] = size(Rotation); 
for i = 1:NX
    disp([num2str(i),' / ',num2str(NX)])
    for j = 1:NY
        cell1 = create_meta2(cell1,layerMap,i,j,Px,Py,a(Img_C(i,j)),b(Img_C(i,j)),Rotation(i,j),lam_SPP(Img_C(i,j)),Phi_grad);
        %创建(i,j)单元的cell
    end
end

topcell = add_ref(topcell,cell1);
gdsl = gds_library('gdsl',topcell,cell1);
write_gds_library(gdsl,'!MetaSurface_rect5.gds');

function [Hol1,Hol2] = holgram_cal()
maxn = 100;% 设置迭代次数
Hol_img1 = imread("C:\Users\De'l'l\Desktop\bnu(1).png");
Hol_img2 = imread("C:\Users\De'l'l\Desktop\bnu(2).png");

Comb_img = zeros([size(Img_color),3]);
Comb_hologram = zeros([size(Img_color)]);

    figure;
    subplot(131)
%     Im = rgb2gray(Hol_img(:,:,:,k));
%     Im = double(Im);
    Im = Hol_img(:,:,k);
    Im = double(Im);
    Im = imresize(Im,size(Img_color),'bilinear');
    im_original = Im./max(max(Im));
    %im_original = 1 - im_original;
    
    imshow(im_original,[]);
    title('理想光场分布');
    
    % im_hologram表示全息面的图像，im_image表示成像面图像
    im_image = im_original.*exp(2i*pi*rand(size(im_original)));% 加入随机相位
    for n = 1:maxn % 迭代到指定次数结束
    % 反向传播计算全息图
    im_hologram = ifft2(ifftshift(im_original.*exp(1i*angle(im_image))));
    % 取全息图相位，振幅全部置一，模拟正向传播，计算像面图像
    im_image = fftshift(fft2(Img(:,:,k).*exp(1i*angle(im_hologram))));
    end
    angle_hologram = Img(:,:,k).*angle(im_hologram);
    angle_hologram = mod(angle_hologram,2*pi);
    % 相位全息图显示
    subplot(132)
    imshow(angle_hologram,[0,2*pi]);
    %img255=uint8((angle_hologram+pi)/2/pi*255);
    %imwrite(img255,'D:\far_smile.png');
    title('器件面相位分布');
    
    % 光场归一化
    im_field = abs(fftshift(fft2(Img(:,:,k).*exp(1i*angle(im_hologram)))));
    im_field = (im_field -min(min(im_field)))./(max(max(im_field)) - min(min(im_field)));
    % 全息像显示
    subplot(133)
    Himg = zeros([size(im_field),3]);
    Himg(:,:,k) = im_field;
    imshow(Himg);
    title('远场光场分布');
    Comb_img(:,:,k) = im_field; %将三通道全息图的还原像组合
    Comb_hologram = Comb_hologram + angle_hologram; %将三通道全息图组合


figure;
imshow(Comb_img);
figure;
Comb_hologram = mod(Comb_hologram,2*pi); %几何相位存储
imshow(Comb_hologram);

end

function [lamSPP] = LS(lam,eps0)  %SPP波长计算，eps0为介质介电常数实数部分
    eps1 = real(LD(lam,'Ag'));     %此函数需要根据材料改变
    lamSPP = lam*sqrt((eps1+eps0)/(eps1*eps0));
end


function [cell1] = create_meta2(cell1,layerMap,index_X,index_Y,Px,Py,a,b,Rot,lam_eff,Phi_grad)
    dy = 0.5; %y方向偏移
    sign = 1;
    if mod((index_Y-1)*Py*Phi_grad,8*pi) < 4*pi
        dx = mod((index_Y-1)*Py*Phi_grad,8*pi)/(2*pi)*lam_eff;
        sign = 1;
    else
        dx = (8*pi - mod((index_Y-1)*Py*Phi_grad,8*pi))/(2*pi)*lam_eff;
        sign = -1;
    end
    Pos0 = [Px*(index_X-1) + dx+ a/2*cosd(Rot) + b/2*sind(Rot), -Py*(index_Y-1)+ a/2*sind(Rot) - b/2*cosd(Rot)];  %单元格左上角点
    Posc(1) = Pos0(1) - a/2*cosd(Rot) - b/2*sind(Rot); %矩形几何中心点横纵坐标计算
    Posc(2) = Pos0(2) - a/2*sind(Rot) + b/2*cosd(Rot);

    Posc = Posc + [sign*lam_eff/2,-dy]; %新矩形的中心坐标  左旋圆偏振光 %换为-lam_eff/2即为右旋圆偏振光，通过改为其他数值可以制备其他形式的偏振光
    Pos1(1) = Posc(1) +a/2*cosd(Rot+90) +b/2*sind(Rot+90); %计算新矩形的固定点
    Pos1(2) = Posc(2) +a/2*sind(Rot+90) -b/2*cosd(Rot+90);

    %Pos1 = Pos0;   %单元格左上角点
    %Pos2 = [Px*index_X, -Py*(index_Y-1)];  %单元格右上角点
    waveguideA = Waveguide(b,layerMap.FullCore,5,3);
    info1 = CursorInfo(Pos0,Rot+180 ,3.17);
    [cell1,~] = PlaceRect(cell1,info1,a,waveguideA.w,waveguideA.layer,waveguideA.dtype);
    %waveguideA = Waveguide(b,layerMap.FullCore,5,3);
    info1 = CursorInfo(Pos1,Rot+270 ,3.17);
    [cell1,~] = PlaceRect(cell1,info1,a,waveguideA.w,waveguideA.layer,waveguideA.dtype);
%     if Rot > 0 && Rot< 90  %锐角旋转的情况（旋转角度不同，几何不同）
%         %Pos1(1) = Pos1(1)+period/2;
% 
%             %平行线计算终止条件
%         while Pos1(1) < Px*index_X || Pos1(2) > -Py*index_Y  
%             Point(1,1:2) = Pos1; %第一个点（与Point3共平行线）
% 
%              %计算下一个平行线的Pos1，分情况讨论
%             if Pos1(1) + period/sind(Rot) <= Pos0(1) + Px  %当下一个Pos1与上一个点在同一个边上时（上边）
%                 Pos1(1) = Pos1(1) + period/sind(Rot);
%             elseif Pos1(1) ~= Pos0(1) + Px %两者不同边，转变情况
%                 Pos1(2) = Pos0(2) - (period - (Pos0(1)+Px-Pos1(1))*sind(Rot))/cosd(Rot);
%                 Pos1(1) = Pos0(1) + Px; %转变后都在竖直边上，所以x坐标为常数
%             else    %转变后，两者在右侧竖直边上
%                 Pos1(2) = Pos1(2) - period/cosd(Rot);
%             end
% 
% 
%             Point(2,1:2) = Pos1; %Point2即为新Pos1坐标点
% 
%             %Point3为过Point1的平行线与边界的另一个交点
%             Point(3,1:2) = [Pos0(1), Point(1,2) - (Point(1,1)-Pos0(1))*tand(Rot)];
%             if Point(3,2) < Pos0(2) - Py
%                 Point(3,1:2) = [Point(1,1)-( Point(1,2)-Pos0(2)+Py )/tand(Rot) , Pos0(2)-Py];
%             end
% 
% %             Point(4,1:2) = [Pos0(1), Pos0(2) - (Point(2,1)-Pos0(1))*tand(Rot)];
% %             if Point(4,2) < Pos0(2) - Py
% %                 Point(4,1:2) = [Point(2,1)-Py/tand(Rot) , Pos0(2)-Py];
% %             end
% 
%             waveguideA = Waveguide(period*Crate,layerMap.FullCore,5,3);
% 
%             info1 = CursorInfo(Point(1,1:2),Rot+180 ,3.17);
%             [cell1,~] = PlaceRect(cell1,info1,sqrt(sum((Point(1,:)-Point(3,:)).^2)),waveguideA.w,waveguideA.layer,waveguideA.dtype);
% 
%         end
%     elseif Rot > 90
%          while Pos2(1) > Pos0(1) || Pos2(2) > -Py*index_Y 
%             Point(1,1:2) = Pos2;
% 
%             if Pos2(1) - period/sind(180-Rot) >= Pos0(1)
%                 Pos2(1) = Pos2(1) - period/sind(180-Rot);
%             elseif Pos2(1) ~= Pos0(1)
%                 Pos2(2) = Pos0(2) - (period - (Pos2(1)-Pos0(1))*sind(180-Rot))/cosd(180-Rot);
%                 Pos2(1) = Pos0(1) ;
%             else
%                 Pos2(2) = Pos2(2) - period/cosd(180-Rot);
%             end
% 
%             Point(2,1:2) = Pos2;
% 
%             Point(3,1:2) = [Pos0(1)+Px , Point(1,2) - (Pos0(1)+Px-Point(1,1)) *tand(180-Rot)];
%             if Point(3,2) < Pos0(2) - Py
%                 Point(3,1:2) = [Point(1,1) + (Point(1,2)-Pos0(2)+Py)/tand(180-Rot) , Pos0(2)-Py];
%             end
% 
% %             Point(4,1:2) = [Pos0(1), Pos0(2) - (Point(2,1)-Pos0(1))*tand(Rot)];
% %             if Point(4,2) < Pos0(2) - Py
% %                 Point(4,1:2) = [Point(2,1)-Py/tand(Rot) , Pos0(2)-Py];
% %             end
% 
%             waveguideA = Waveguide(period*Crate,layerMap.FullCore,5,3);
% 
%             info1 = CursorInfo(Point(1,1:2),Rot+180,3.17);
%             [cell1,~] = PlaceRect(cell1,info1,sqrt(sum((Point(1,:)-Point(3,:)).^2)),waveguideA.w,waveguideA.layer,waveguideA.dtype);
% 
%          end

%    end


end
    