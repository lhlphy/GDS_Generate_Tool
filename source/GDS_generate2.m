%Versatile full-colour nanopainting enabled 
%by a pixelated plasmonic metasurface
%Code Author:Li Haolin
%Last revised date:2024/11/5
%矩形块阵列


% set the work path to parentFolder
currentFile = mfilename('fullpath');  % 获取当前脚本或文件的完整路径
currentFolder = fileparts(currentFile); % 获取文件所在文件夹的路径
parentFolder = fileparts(currentFolder);  % 获取上一层文件夹的路径
cd([parentFolder,'\Generated']); % 设置工作目录为上一层文件夹

[log,cad,cellname,topcell,layerMap] = InitializeCell();
refs = [];
refs = GetRefsFloorplan(refs);

Px = 0.550; %大单元格x周期
Py = 0.550; %大单元格y周期
a = 0.2*[0,2.35,1.85,1.40]; %矩形块长轴
b = 0.06*[0,2.35,1.85,1.40]; %矩形块短轴
Rotation = Comb_hologram*180/(2*pi); %注意几何相位为转动角Rot的两倍
imgmax = max(max(Img_color)); %色彩匹配编号
Img_color = Img_color +(4-imgmax);  %将0，1，2，3转为1，2，3，4

cell1 = gds_structure('cell1');
%开始逐个单元绘制版图
[NX,NY] = size(Rotation); 
for i = 1:NX
    disp([num2str(i),' / ',num2str(NX)])
    for j = 1:NY
        cell1 = create_meta2(cell1,layerMap,i,j,Px,Py,a(Img_color(i,j)),b(Img_color(i,j)),Rotation(i,j));
        %创建(i,j)单元的cell
    end
end

topcell = add_ref(topcell,cell1);
gdsl = gds_library('gdsl',topcell,cell1);
write_gds_library(gdsl,'!MetaSurface_rect2.gds');


function [cell1] = create_meta2(cell1,layerMap,index_X,index_Y,Px,Py,a,b,Rot)
    Pos0 = [Px*(index_X-1)+ a/2*cosd(Rot) + b/2*sind(Rot), -Py*(index_Y-1)+ a/2*sind(Rot) - b/2*cosd(Rot)];  %中心坐标
    %Pos1 = Pos0;   %单元格左上角点
    %Pos2 = [Px*index_X, -Py*(index_Y-1)];  %单元格右上角点
    waveguideA = Waveguide(b,layerMap.FullCore,5,3); %生成纳米槽（矩形）单元
    info1 = CursorInfo(Pos0,Rot+180 ,3.17); %提供单元在cell上的位置信息
    [cell1,~] = PlaceRect(cell1,info1,a,waveguideA.w,waveguideA.layer,waveguideA.dtype); %放置纳米槽到cell上
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
    