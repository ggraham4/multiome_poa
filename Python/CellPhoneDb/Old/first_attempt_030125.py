import anndata 
import seaborn  
import pandas as pd

#I am going to kill myself
import cellphonedb

import anndata as ad
###OH MY FUCKING GOD IT FINALLY WORKED THANK YOU FUCKING JESUS IN MY ASS IDK WHAT I DID I THINK I HAD TO SPECIFY THE PATH WHERE PIP INSTALLED IT BUT THAT STILL DIDNT WORK
### IN ANY CASE, ALWAYS USE PYTHON 3.12 IDGAF WHAT YOU HAVE TO DO

##the command is pip install [package] --target /Users/ggraham/Library/Python/3.11/lib/python/site-packages --upgrade

###ok I need to make an anndata object I know that 

obj = anndata.io.read_h5ad('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA_object_anndata.h5ad')
###IT WORKED

import sysconfig
print(sysconfig.get_paths()["purelib"])