## Before using...
Rely on MatlabGDSPhotonicsToolbox, install `MatlabGDSPhotonicsToolbox` first before using.  
It is recommended to use the third-party toolkit installation tool that comes with MATLAB to install it, and search for "GDS" to find `MatlabGDSPhotonicsToolbox`  
## Userguide
This project is designed to generate GDS layouts of metasurfaces that simultaneously achieve full-color display and multi-channel holography (RGB three channels), for use by researchers in micro-nano fabrication. The project takes in two images, one for full-color display and one for multi-channel holography. The design for full-color display is based on the matching of structure and resonant wavelength (using a specific material, Ag), while the multi-channel holography is realized through the GS algorithm, utilizing the geometric phase of the structure.  
#### The structures provided by this project include:  
*``Surface shallow grating arrays  
Nano-slot array  
Orthogonal double nano-slot array  
Orthogonal double nano-slot array with phase modulation propagation (uniform phase gradient)  
Dual-zone metasurface  
Three-zone metasurface``*
#### Important Functionality Description  
1. In this tool, we use `metasurface_calculate.m` to transform pictures to info matrix (shape, rotation info for .gds).Then using the calculated info to generate .gds files.
**So, run `metasurface_calculate.m` first, then run `GDS_generate<i>.m` to generate .gds files**
2. You can modify the files or paths of `Origin_picture\pic.png` and `Origin_picture\pic(1).png` (in `metasurface_calculate.m`) to change the images for full-color display and multi-channel holographic display, respectively.
3. You can change the image size in metasurface_calculate.m by modifying `imgsize`  
`imgsize = [100,100]; % 图片大小（像素单元数量），过大的图片会导致生成过慢`  
A larger image will cause slower generation speed.
4. If you wish to use different materials, please modify the following parameters in `GDS_generate<i>.m`:
```
Px = 0.550; %大单元格x周期
Py = 0.550; %大单元格y周期
a = 0.2*[0,2.35,1.85,1.40]; %矩形块长轴
b = 0.06*[0,2.35,1.85,1.40]; %矩形块短轴
```
In `GDS_generate1.m`, modify the following parameters:
```
Px = 2; %大单元格x周期 <2.6
Py = 2; %大单元格y周期 <2.6
period = [0,0.58,0.4833,0.4143]; %单元格内光栅周期 
Crate = [0.5,380/580,240/483.3,185/414.3]; %光栅占空比
```
The last three parameters in the list correspond to the structural parameters for when resonance occurs at the three different wavelengths for RGB (displaying different colors). For any material, these parameters can be determined using one of the following three methods:  
- Consult literature: If there are existing metasurface designs using the same material, use the same parameters.  
- Use Comsol simulation: Perform a parametric scan of different structural parameters to find the optimal parameters corresponding to the required resonant wavelengths.  
- Use the dielectric model of the material for theoretical calculations: Obtain the analytical relationship between the structural parameters and resonant wavelengths (if you can...).

5. If you do not wish to use images to generate GDS information, you can also customize the GDS information matrix. Please directly define the `Img_color` matrix (which determines the structural parameters and colors) and the `Rotation` matrix (which determines the rotation angles, geometric phase, and is directly related to the `Comb_hologram` matrix) in `GDS_generate<i>.m`. It is recommended to define both matrices before running `GDS_generate<i>.m`. Please ensure that the sizes of both matrices match.  
**if you decide to define `Rotation` matrix and `Img_color` matrix by yourself, please \*delete\* the following line:**  
`Rotation = Comb_hologram*180/(2*pi); %注意几何相位为转动角Rot的两倍`
## File List:  
#### Folder: source  
- metasurface_calculate.m: Used to calculate the geometric parameters of the metasurface, which are provided to GDS_generate.m for layout generation.  
- GDS_generate1.m: Uses geometric parameters to generate EBL layout (GDS file) for drawing surface shallow grating arrays.  
- GDS_generate2.m: Generates a nano-slot array layout.  
- GDS_generate3.m: Generates an orthogonal double nano-slot array layout.  
- GDS_generate4.m: Generates an orthogonal double nano-slot array layout with phase modulation propagation (uniform phase gradient).  
- GDS_generate5.m: Generates a dual-zone metasurface layout used for separating left and right circularly polarized light.  
- GDS_generate6.m: Generates a three-zone metasurface layout used for multi-holography (3x3=9) display.  
- LD.m: A function used in GDS_generate4-6 to calculate SPP wavelengths, with a built-in Drude-Lorentz material library (only includes several commonly used materials).  

#### Folder: Original_picture  
- Original_image.png: Original image for full-color display and multi-channel holography.  
- pic.png: Preset original image for full-color display.  
- pic(1).png: Preset original image for multi-channel holography.  

#### Folder: Generated  
- Full_Color_Display.png: Full-color display image after metasurface color matching.  
- metasurface_colourful_cells.png: Color map of metasurface structural units, where each pixel corresponds to a metasurface unit color.  
- MetaSurface_rect<i>: GDS files (metasurface layouts) output from GDS_generate<i>.m, i=1,2,3,4,5,6. After fine-tuning (dose adjustment), these files are ready for EBL fabrication.  

### 文件列表：  
#### Folder: source  
- metasurface_calculate.m ：用于计算超表面的几何参数，提供给GDS_generate.m绘制版图  
- GDS_generate1.m ：利用几何参数绘制EBL版图（GDS文件），绘制表面浅光栅阵列  
- GDS_generate2.m ：绘制纳米槽阵列  
- GDS_generate3.m ：绘制正交双纳米槽阵列  
- GDS_generate4.m ：绘制传播相位调制的正交双纳米槽阵列（均匀相位梯度）  
- GDS_generate5.m ：绘制双分区超表面，用于将左右旋偏振光分离  
- GDS_generate6.m ：绘制三分区超表面，用于多全息图（3*3=9）显示
- LD.m：GDS_generate4-6中需要的用于计算SPPs波长的函数，内置Drude-Lorentz模型材料库（仅包含几种常用的材料）  
#### Folder: Original_picture  
- Original_image.png ：用于全彩显示和多通道全息的原始图像  
- pic.png : 预设的用于全彩显示的原始图像  
- pic(1).png : 预设的用于多通道全息的原始图像  
#### Folder: Generated  
- 全彩显示图.png ：超表面色彩匹配后的全彩显示图  
- metasurface_colourful_cells.png ：超表面结构单元色彩图，每个像素点对应一个超表面单元色彩  
- MetaSurface_rect<i>：GDS_generate<i>.m，i=1,2,3,4,5,6运行后输出的GDS文件（超表面版图），经过微调（剂量调整）即可用于EBL加工  
