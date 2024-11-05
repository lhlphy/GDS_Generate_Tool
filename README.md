### File List:  
#### Folder: source  
- metasurface_calculate.m: Used to calculate the geometric parameters of the metasurface, which are provided to GDS_generate.m for layout generation.  
- GDS_generate.m: Uses geometric parameters to generate EBL layout (GDS file) for drawing surface shallow grating arrays.  
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
- GDS_generate.m ：利用几何参数绘制EBL版图（GDS文件），绘制表面浅光栅阵列  
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
