#  hiviewer （old name: Fits_Inspect_Overlap)

This code is used to: Inspect and show the fits image,  Overlap different images (or contours), Process the data cube preliminarily.

Please read README and Fits_example.ipynb in doc directory.

## Notes

* 此射电图像包的一个简要介绍在[这里](https://zhuanlan.zhihu.com/p/595278094).

* 欢迎关注astroR2的[zhihu](https://www.zhihu.com/people/stellarxu).

## Logs

* Newest version: v1.3 on 2023/01/01
* Updated on 2023/01/01
    * 增加pv图。修改了一些bug

* Updated on 2022/09/03
    * 更改包的名字为hiviewer,这个过程可能产生bug

* Updated on 2021/01/26
    * 修改了一点默认参数，谱线可以按pixel抽取
    * 修改了把nan补0的bug
    * moment 2 支持三种type

            'variance': M2

            'sigma'：sigma = np.sqrt(M2)
            
            'fwhm':fwhm = np.sqrt(8ln2)*sigma
    * tag v1.2

* Updated on 2021/01/10
    * 增加了三种reproject方法
        * [Montage](https://montage-wrapper.readthedocs.io/en/latest/)
        * [Kapteyn](https://www.astro.rug.nl/software/kapteyn/maputilstutorial.html#re-projections-and-image-overlays)
        * [reproject](https://reproject.readthedocs.io/en/stable/celestial.html#adaptive-resampling)
    * 增加了椭圆拟合一个星系，画沿长轴的cube截面
    * 更改了moudule结构，更合理
    
* ** IMPORTANT WARNING **   on 2020/12/11
    
    如果使用了reproject，重新投影，改变pixel的大小、分辨率，单位一定要慎重！

    确保单位是**surface brightness**，a unit of flux per area，比如Jy/beam,Jy/sr,W/m^2/sr等，而不是flux，比如mJy,Jy/pixel,W/m^2等。如果是flux，you will need to scale by a factor of  (the new pixel area)/(the old pixel area)。

    具体说明：http://research.endlessfernweh.com/convolution-and-regridding/  

* [Kapteyn Package for Python 3](https://www.astro.rug.nl/software/kapteyn/index.html) 是一个功能更加齐全的包，download the file [kapteyn-3.0.tar.gz](https://www.astro.rug.nl/software/kapteyn/kapteyn-3.0.tar.gz)。

    v1.1之后的版本我会把Kapteyn、[reproject](https://reproject.readthedocs.io/en/stable/celestial.html#adaptive-resampling)中的convolve_reproject功能替换出来，以免Montage的方法不好安装。他们同样不适用于flux。 on 2020/12/1

