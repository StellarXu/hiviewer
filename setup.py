import setuptools
import glob
# with open("README.md", "r") as fh:
#     long_description = fh.read()


setuptools.setup(
    name="hiviewer", # Replace with your own username
    version="1.3.1	",
    author="Xu Chen etc.",
    author_email="xuchen@nao.cas.cn",
    description="A Python module to: Check and show the fits image,  Overlap different images (or contours), Process the data cube preliminarily.",
#     long_description=long_description,
#     long_description_content_type="text/markdown",
#     url="",
    packages=['hiviewer'],
    
    install_requires=['numpy>=1.12','matplotlib','scipy','spectral-cube>=0.6.0','astropy'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Linux",
    ],
    python_requires='>=3.6',
    zip_safe=False,
)
