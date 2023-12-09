README v1.3

This code is used to: Inspect and show the fits image,  Overlap different images (or contours), Process the data cube preliminarily.

# 关于 hiviewer 模块
* 用于
  * 使用Python为cube出图，调出Fits文件的信息，叠加不同fits的图

* 包含
  * core.py                        主要的class和操作  (./hiviewer/)
  * overlap.py                     叠图的函数 (./hiviewer/)
  * HIsource_fio.py                椭圆拟合一个HI源的HIviewer简化版  (./hiviewer/)
  * auto_compare.py                一个可以通过终端自动画出图片的程序  (./hiviewer/)
  * compare_instance.py            与上一个对应的半自动版本程序  (./hiviewer/)
  * utils.py                       一些函数 (./hiviewer/)
  * Fits_example.ipynb             非常详细的notebook实例  (./doc/)
  * MESSIER_033_I_103aE_dss1.fits  Palomar天文台拍摄的光学波段的M33  (./data/)

* 需要的Python包：
    * Python3环境
    * numpy,matplotlib,astropy,spectral-cube
    * 重新投影需要montage_wrapper(它又需要sudo apt get montage)或者kapteyn或者reproject包
    * 寻找极值点位置需要skimage

* 如果有问题或者bug如何联系我？
  xuchen@nao.cas.cn

# 安装
使用setup.py

  下载安装包后 ** 先 cd 切换到代码(setup.py)所在目录下 **
  * 安装
    ```
    python setup.py install -f
    ```
  * 卸载
    ```
    python -c 'import hiviewer;print(hiviewer.__file__)'
    ```
    找到安装后的文件夹，删除。

# 示例

## 准备

从这里下载其他示例data
* FAST对M33 HI的一个数据cube，[点击这里下载](https://www.aliyundrive.com/s/Xs5M5moUwjp),提取码：34gm
* Arecibo data: <http://www.naic.edu/~ages/public_cubes/M33_local.fits.gz>，article: <http://adsabs.harvard.edu/abs/2016MNRAS.456..951K>
* VLA data and article: <https://www.researchgate.net/publication/260025915_M33_HI_maps>

1. 在任意位置新建文件夹 test 
2. 复制doc里的ipynb 到 test 文件夹下 
3. 在test下新建文件夹 data, fig，用来存放数据和图片，示例数据放在data文件夹下
4. 学习Fits_example.ipynb 即可。使用jupyter notebook再配合DS9/CARTA使用更佳。

## 关于class FitsPic 
这个类的属性与方法
```
    Attributes & Properties
    ----------
    filename: str   Name of FITS file
    ndim:  dimension number of the data
    ctype: ['coor--pro']    Coordinate type  and  projection.
    crval: Coordinate value at reference point (one per axis).
    crpix: Array location of the reference point in pixels (one per axis)
    cuint: Units for stated coordinate axes
    cdelt: Coordinate increment at reference point
    velo:  spectral axis of the cube 
    fill_nans_data: fill Nans of the data with zero.
    beam:  resolution information of the telescope
    
    Methods
    ----------
    deg2pix:         ra,dec coordinates to pixel position 
    pix2deg:         pixel position to ra,dec coordinates
    velo2sel:        velocity(km/s) to selected slice number
    get/plot_slice():get a 2D slice of the cube / show the  slice image.
    get/plot_spec(): get/plot the spectrum of a certain point.
    plot_side_slice: show a cube slice from one side.
    slab_moment():   compute 3 moment maps.
    plot_contour:    plot the contour on a certain slice
    reproj_convolve_2D: astropy.convolve and reproject the 2D fits to target fits
    reproj_convolve_cube: SpectralCube.convolve and reproject the fits cube to target fits.
    center_levels :  Expirical contour levels for a 2D FitsPic central part.
    
```
以上以及overlap操作的具体说明请直接help（  ）。



## Fits_example.ipynb             
这是一个非常详细的notebook实例，你可以了解FitsPic类的各种属性和方法，以及如何使用此包。

主要功能有：
* 演示FitsPic类的基础属性和方法
* 制作moment maps
* 绘制截面，pv图，等高线图，查看谱线等
* 卷积与重新投影

```
from hiviewer import FitsPic
from hiviewer.overlap import overlap_img_contour,overlap_two_contours,two_contours_with_peak
```
如果按照前面的安装方法，是可以按照上面方法import的。如果安装不成功或者找不到路径，可以在同一目录下，临时添加路径
```
import sys
sys.path.append('../../../pkgs/hiviewer-master')#这是包的路径，例如前面的/test/路径
from hiviewer import *
```
脚本同理

画图的话，我强烈推荐jupyter notebook,除非你的服务器上没有，那么可以考虑终端内运行python脚本。因为core.py这个程序里有许多参数是可调的，有一些常用参数可以外部更改，其他的可以直接在下载路径里修改core.py这个脚本。例如
```
 fig = plt.figure(figsize=(10,10))#调整大小
 ax.set_xlabel('设置标签',fontsize=16)#标签字号
 ax.set_xlim(1400，1405)#设置范围
 ax.imshow(data,vmin=vmin_max[0],vmax=vmin_max[1],alpha=0.5,cmap='rainbow')#上下限、透明度、colormap
 ax.tick_params(labelsize=17) #轴上数字字号
```
等等，诸如此类的。

具体说明见notebook。



# 如何找到多波段数据？
1. 推荐先到[NED，NASA Extragalactic Database](https://ned.ipac.caltech.edu/)上搜索一下,找一下Image里有哪些。
[images search](http://ned.ipac.caltech.edu/forms/images.html)
[data search](http://ned.ipac.caltech.edu/forms/data.html>)
在后面的两个链接里，可以在最后找到EXTERNAL ARCHIVES AND SERVICES for MESSIER 031 Help一栏，也会给出很多获得数据的sites，resources的连接，或者直接retrieve到.

当然也推荐[CDS](https://cds.u-strasbg.fr/)数据库里搜索也很好，例如寻找EBHIS的一些[数据](http://cdsarc.u-strasbg.fr/viz-bin/qcat?J/A+A/585/A41)

2. [Lambda,NASA](https://lambda.gsfc.nasa.gov/product/)给出了很多data product
比如[HI Surveys、连续谱的、偏振的](https://lambda.gsfc.nasa.gov/product/foreground/fg_diffuse.cfm)从中可以找到某些巡天、论文或者data下载路径，虽然他叫Legacy Archive for Microwave Background Data Analysis...

[heasarc,NASA](https://heasarc.gsfc.nasa.gov/)上也有很多，比如首页就有查询xamin

3. 特意寻找一些巡天，常用的比如
[ALFALFA](http://egg.astro.cornell.edu/alfalfa/data/index.php)Arecibo
[Arecibo Galaxy Environment Survey,AGES](http://www.naic.edu/~ages/)
[The VLA FIRST Survey](http://sundog.stsci.edu/)中性氢巡天
[The NRAO VLA Sky Survey ,NVSS](https://www.cv.nrao.edu/nvss/)他在首页里也提到了很多Other large-scale radio surveys which may be of interest include
[The Effelsberg-Bonn HI Survey (EBHIS)](https://astro.uni-bonn.de/~jkerp/index.php)或者[这个马普所网页](https://www.mpifr-bonn.mpg.de/pressreleases/2015/9)
[GALFA HI Data](https://purcell.ssl.berkeley.edu/)
帕克斯中性氢巡天[HIPASS](http://hipass.anu.edu.au/)与焦德雷班克中性氢巡天HIJASS等等
```
Table B.1 Comparison of major blind HI surveys
Survey	Area	Beam	Vmax	Vresa	ts	rmsb	Ndet	min MHI c	Ref
 	(deg2)	(arcmin)	(km/s)	(km/s)	(s)	(mJy)	 	(Msun)	 
AHISS	65	 3.3	-700 - 7400	16	var	0.7	65	1.9x106	1
ADBS	430	 3.3	-650 - 7980	34	12	3.6	265	9.9x106	2
WSRT	1800	49. 	-1000 - 6500	17	60	18	155	4.9x107	3
Nancay CVn	800	4 x 20	-350 - 2350	10	80	7.5	33	2.0x107	4
HIJASS	1115	12. 	-1000 - 10000d	18	400	13	222	3.6x107	5
HIJASS-VIR	32	12. 	500 - 2500	18	3500	4.	31	1.1x107	6
HIDEEP	60	15.5	-1280 - 12700	18	9000	3.2	173	8.8x106	7
HIZSS	1840	15.5	-1280 - 12700	27	200	15.	110	4.1x107	8
HICAT	21341	15.5	300 - 12700	18	450	13.	4315	3.6x107	9
HIPASS	 	15.5	300 - 12700	18	450	13.	(6000)	3.6x107	10
AUDS	0.4	 3.5	-960 - 47000e	TBD	70 × 3600	0.02	(40)	0.6x10 6	11
AGES	TBD	 3.5	-960 - 47000e	TBD	300	0.5	TBD	1.4x106	12
ALFALFA	7000	 3.5	-2000 - 18000	11	28	1.6	(30000)	4.4x106	
```

只是一部分，不全

4. 去搜索论文，看看有没有告诉你数据在哪。
比如我在ReseaechGate上找到的VLA的M33


end
2020/06/02
2020/12/03
2023/01/01


