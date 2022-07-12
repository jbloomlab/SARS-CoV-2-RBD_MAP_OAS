# Aggregate variant counts for all samples
Separate `Snakemake` rules count the observations of each variant in each sample from the Illumina barcode sequencing.
This Python Jupyter notebook aggregates all of this counts, and then adds them to a codon variant table.

## Set up analysis
### Import Python modules.
Use [plotnine](https://plotnine.readthedocs.io/en/stable/) for ggplot2-like plotting.

The analysis relies heavily on the Bloom lab's [dms_variants](https://jbloomlab.github.io/dms_variants) package:


```python
import glob
import itertools
import math
import os
import warnings

import Bio.SeqIO

import dms_variants.codonvarianttable
from dms_variants.constants import CBPALETTE
import dms_variants.utils
import dms_variants.plotnine_themes

from IPython.display import display, HTML

import pandas as pd

from plotnine import *

import yaml
```

Set [plotnine](https://plotnine.readthedocs.io/en/stable/) theme to the gray-grid one defined in `dms_variants`:


```python
theme_set(dms_variants.plotnine_themes.theme_graygrid())
```

Versions of key software:


```python
print(f"Using dms_variants version {dms_variants.__version__}")
```

    Using dms_variants version 1.4.0


Ignore warnings that clutter output:


```python
warnings.simplefilter('ignore')
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Make output directory if needed:


```python
os.makedirs(config['counts_dir'], exist_ok=True)
```

## Read barcode counts / fates
Read data frame with list of all samples (barcode runs):


```python
print(f"Reading list of barcode runs from {config['barcode_runs']}")

barcode_runs = (pd.read_csv(config['barcode_runs'])
                .assign(sample_lib=lambda x: x['sample'] + '_' + x['library'],
                        counts_file=lambda x: config['counts_dir'] + '/' + x['sample_lib'] + '_counts.csv',
                        fates_file=lambda x: config['counts_dir'] + '/' + x['sample_lib'] + '_fates.csv',
                        )
                .drop(columns='R1')  # don't need this column, and very large
                )

assert all(map(os.path.isfile, barcode_runs['counts_file'])), 'missing some counts files'
assert all(map(os.path.isfile, barcode_runs['fates_file'])), 'missing some fates files'

display(HTML(barcode_runs.to_html(index=False)))
```

    Reading list of barcode runs from data/barcode_runs.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>date</th>
      <th>experiment</th>
      <th>target</th>
      <th>secondary</th>
      <th>library</th>
      <th>antibody</th>
      <th>concentration</th>
      <th>sort_bin</th>
      <th>selection</th>
      <th>sample</th>
      <th>experiment_type</th>
      <th>number_cells</th>
      <th>frac_escape</th>
      <th>HutchBase_name</th>
      <th>sample_lib</th>
      <th>counts_file</th>
      <th>fates_file</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>220623</td>
      <td>VicOAS_1-14</td>
      <td>Wuhan-Hu-1</td>
      <td>none</td>
      <td>lib1</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>reference</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>ab_selection</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>220623_Wuhan_lib1_ref</td>
      <td>none-Wuhan-Hu-1-none-ref_lib1</td>
      <td>results/counts/none-Wuhan-Hu-1-none-ref_lib1_counts.csv</td>
      <td>results/counts/none-Wuhan-Hu-1-none-ref_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_1-14</td>
      <td>Wuhan-Hu-1</td>
      <td>none</td>
      <td>lib2</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>reference</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>ab_selection</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>220623_Wuhan_lib2_ref</td>
      <td>none-Wuhan-Hu-1-none-ref_lib2</td>
      <td>results/counts/none-Wuhan-Hu-1-none-ref_lib2_counts.csv</td>
      <td>results/counts/none-Wuhan-Hu-1-none-ref_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_1-14</td>
      <td>BA.1</td>
      <td>none</td>
      <td>lib1</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>reference</td>
      <td>none-BA.1-none-ref</td>
      <td>ab_selection</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>220623_BA1_lib1_ref</td>
      <td>none-BA.1-none-ref_lib1</td>
      <td>results/counts/none-BA.1-none-ref_lib1_counts.csv</td>
      <td>results/counts/none-BA.1-none-ref_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_1-14</td>
      <td>BA.1</td>
      <td>none</td>
      <td>lib2</td>
      <td>none</td>
      <td>0</td>
      <td>ref</td>
      <td>reference</td>
      <td>none-BA.1-none-ref</td>
      <td>ab_selection</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>220623_BA1_lib2_ref</td>
      <td>none-BA.1-none-ref_lib2</td>
      <td>results/counts/none-BA.1-none-ref_lib2_counts.csv</td>
      <td>results/counts/none-BA.1-none-ref_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_01</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib1</td>
      <td>WOO-3</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1315751.0</td>
      <td>0.328</td>
      <td>VicOAS01_lib1_abneg</td>
      <td>WOO-3-BA.1-FLAG-abneg_lib1</td>
      <td>results/counts/WOO-3-BA.1-FLAG-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-3-BA.1-FLAG-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_01</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib2</td>
      <td>WOO-3</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1420019.0</td>
      <td>0.345</td>
      <td>VicOAS01_lib2_abneg</td>
      <td>WOO-3-BA.1-FLAG-abneg_lib2</td>
      <td>results/counts/WOO-3-BA.1-FLAG-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-3-BA.1-FLAG-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_02</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib1</td>
      <td>WOO-4</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>911578.0</td>
      <td>0.226</td>
      <td>VicOAS02_lib1_abneg</td>
      <td>WOO-4-BA.1-FLAG-abneg_lib1</td>
      <td>results/counts/WOO-4-BA.1-FLAG-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-4-BA.1-FLAG-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_02</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib2</td>
      <td>WOO-4</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>853762.0</td>
      <td>0.226</td>
      <td>VicOAS02_lib2_abneg</td>
      <td>WOO-4-BA.1-FLAG-abneg_lib2</td>
      <td>results/counts/WOO-4-BA.1-FLAG-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-4-BA.1-FLAG-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_03</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib1</td>
      <td>WOO-5</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1092980.0</td>
      <td>0.299</td>
      <td>VicOAS03_lib1_abneg</td>
      <td>WOO-5-BA.1-FLAG-abneg_lib1</td>
      <td>results/counts/WOO-5-BA.1-FLAG-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-5-BA.1-FLAG-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_03</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib2</td>
      <td>WOO-5</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1206355.0</td>
      <td>0.299</td>
      <td>VicOAS03_lib2_abneg</td>
      <td>WOO-5-BA.1-FLAG-abneg_lib2</td>
      <td>results/counts/WOO-5-BA.1-FLAG-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-5-BA.1-FLAG-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_04</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib1</td>
      <td>WOO-6</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1153142.0</td>
      <td>0.285</td>
      <td>VicOAS04_lib1_abneg</td>
      <td>WOO-6-BA.1-FLAG-abneg_lib1</td>
      <td>results/counts/WOO-6-BA.1-FLAG-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-6-BA.1-FLAG-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_04</td>
      <td>BA.1</td>
      <td>FLAG</td>
      <td>lib2</td>
      <td>WOO-6</td>
      <td>50</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-BA.1-FLAG-abneg</td>
      <td>ab_selection</td>
      <td>1197316.0</td>
      <td>0.294</td>
      <td>VicOAS04_lib2_abneg</td>
      <td>WOO-6-BA.1-FLAG-abneg_lib2</td>
      <td>results/counts/WOO-6-BA.1-FLAG-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-6-BA.1-FLAG-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_05</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-3</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1618434.0</td>
      <td>0.384</td>
      <td>VicOAS05_lib1_abneg</td>
      <td>WOO-3-BA.1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-3-BA.1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-3-BA.1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_05</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-3</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1526528.0</td>
      <td>0.400</td>
      <td>VicOAS05_lib2_abneg</td>
      <td>WOO-3-BA.1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-3-BA.1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-3-BA.1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_06</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-4</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1440151.0</td>
      <td>0.351</td>
      <td>VicOAS06_lib1_abneg</td>
      <td>WOO-4-BA.1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-4-BA.1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-4-BA.1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_06</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-4</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1439159.0</td>
      <td>0.371</td>
      <td>VicOAS06_lib2_abneg</td>
      <td>WOO-4-BA.1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-4-BA.1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-4-BA.1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_07</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-5</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1347197.0</td>
      <td>0.355</td>
      <td>VicOAS07_lib1_abneg</td>
      <td>WOO-5-BA.1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-5-BA.1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-5-BA.1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_07</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-5</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1402918.0</td>
      <td>0.369</td>
      <td>VicOAS07_lib2_abneg</td>
      <td>WOO-5-BA.1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-5-BA.1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-5-BA.1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_08</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-6</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1465522.0</td>
      <td>0.378</td>
      <td>VicOAS08_lib1_abneg</td>
      <td>WOO-6-BA.1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-6-BA.1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-6-BA.1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_08</td>
      <td>BA.1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-6</td>
      <td>200</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-BA.1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>1624568.0</td>
      <td>0.396</td>
      <td>VicOAS08_lib2_abneg</td>
      <td>WOO-6-BA.1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-6-BA.1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-6-BA.1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_09</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WWW-1</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WWW-1-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>855573.0</td>
      <td>0.204</td>
      <td>VicOAS09_lib1_abneg</td>
      <td>WWW-1-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WWW-1-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WWW-1-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_09</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WWW-1</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WWW-1-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>941787.0</td>
      <td>0.224</td>
      <td>VicOAS09_lib2_abneg</td>
      <td>WWW-1-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WWW-1-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WWW-1-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_10</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WWW-2</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WWW-2-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>798404.0</td>
      <td>0.183</td>
      <td>VicOAS10_lib1_abneg</td>
      <td>WWW-2-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WWW-2-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WWW-2-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_10</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WWW-2</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WWW-2-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>821424.0</td>
      <td>0.191</td>
      <td>VicOAS10_lib2_abneg</td>
      <td>WWW-2-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WWW-2-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WWW-2-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_11</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-3</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>673279.0</td>
      <td>0.164</td>
      <td>VicOAS11_lib1_abneg</td>
      <td>WOO-3-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-3-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-3-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_11</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-3</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-3-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>740321.0</td>
      <td>0.172</td>
      <td>VicOAS11_lib2_abneg</td>
      <td>WOO-3-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-3-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-3-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_12</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-4</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>841108.0</td>
      <td>0.219</td>
      <td>VicOAS12_lib1_abneg</td>
      <td>WOO-4-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-4-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-4-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_12</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-4</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-4-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>960379.0</td>
      <td>0.230</td>
      <td>VicOAS12_lib2_abneg</td>
      <td>WOO-4-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-4-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-4-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_13</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-5</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>697012.0</td>
      <td>0.158</td>
      <td>VicOAS13_lib1_abneg</td>
      <td>WOO-5-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-5-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-5-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_13</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-5</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-5-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>736009.0</td>
      <td>0.171</td>
      <td>VicOAS13_lib2_abneg</td>
      <td>WOO-5-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-5-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-5-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_14</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib1</td>
      <td>WOO-6</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>741465.0</td>
      <td>0.177</td>
      <td>VicOAS14_lib1_abneg</td>
      <td>WOO-6-Wuhan-Hu-1-Strep-abneg_lib1</td>
      <td>results/counts/WOO-6-Wuhan-Hu-1-Strep-abneg_lib1_counts.csv</td>
      <td>results/counts/WOO-6-Wuhan-Hu-1-Strep-abneg_lib1_fates.csv</td>
    </tr>
    <tr>
      <td>220623</td>
      <td>VicOAS_14</td>
      <td>Wuhan-Hu-1</td>
      <td>Strep</td>
      <td>lib2</td>
      <td>WOO-6</td>
      <td>1000</td>
      <td>abneg</td>
      <td>escape</td>
      <td>WOO-6-Wuhan-Hu-1-Strep-abneg</td>
      <td>ab_selection</td>
      <td>782764.0</td>
      <td>0.190</td>
      <td>VicOAS14_lib2_abneg</td>
      <td>WOO-6-Wuhan-Hu-1-Strep-abneg_lib2</td>
      <td>results/counts/WOO-6-Wuhan-Hu-1-Strep-abneg_lib2_counts.csv</td>
      <td>results/counts/WOO-6-Wuhan-Hu-1-Strep-abneg_lib2_fates.csv</td>
    </tr>
  </tbody>
</table>


Confirm sample / library combinations unique:


```python
assert len(barcode_runs) == len(barcode_runs.groupby(['sample', 'library']))
```

Make sure the the libraries for which we have barcode runs are all in our variant table:


```python
# unknown_libs = set(barcode_runs['library']) - set(variants.libraries)
# if unknown_libs:
#     raise ValueError(f"Libraries with barcode runs not in variant table: {unknown_libs}")
```

Now concatenate the barcode counts and fates for each sample:


```python
counts = pd.concat([pd.read_csv(f) for f in barcode_runs['counts_file']],
                   sort=False,
                   ignore_index=True)

print('First few lines of counts data frame:')
display(HTML(counts.head().to_html(index=False)))

fates = pd.concat([pd.read_csv(f) for f in barcode_runs['fates_file']],
                  sort=False,
                  ignore_index=True)

print('First few lines of fates data frame:')
display(HTML(fates.head().to_html(index=False)))
```

    First few lines of counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>barcode</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>GAAACAAATTTCTATA</td>
      <td>5871</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>GACTATCGAATTATTG</td>
      <td>5120</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>TTAAAAGCAACCGACC</td>
      <td>5078</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>CGTAACATTTACATAT</td>
      <td>4863</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>AATGAACCCAAAACAA</td>
      <td>4775</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
  </tbody>
</table>


    First few lines of fates data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>fate</th>
      <th>count</th>
      <th>library</th>
      <th>sample</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>valid barcode</td>
      <td>24796153</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>low quality barcode</td>
      <td>12230713</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>invalid barcode</td>
      <td>7411766</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>unparseable barcode</td>
      <td>817674</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
    <tr>
      <td>failed chastity filter</td>
      <td>0</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
    </tr>
  </tbody>
</table>


## Examine fates of parsed barcodes
First, we'll analyze the "fates" of the parsed barcodes.
These fates represent what happened to each Illumina read we parsed:
 - Did the barcode read fail the Illumina chastity filter?
 - Was the barcode *unparseable* (i.e., the read didn't appear to be a valid barcode based on flanking regions)?
 - Was the barcode sequence too *low quality* based on the Illumina quality scores?
 - Was the barcode parseable but *invalid* (i.e., not in our list of variant-associated barcodes in the codon variant table)?
 - Was the barcode *valid*, and so will be added to variant counts.
 
First, we just write a CSV file with all the barcode fates:


```python
fatesfile = os.path.join(config['counts_dir'], 'barcode_fates.csv')
print(f"Writing barcode fates to {fatesfile}")
fates.to_csv(fatesfile, index=False)
```

    Writing barcode fates to results/counts/barcode_fates.csv


Next, we tabulate the barcode fates in wide format:


```python
display(HTML(fates
             .pivot_table(columns='fate',
                          values='count',
                          index=['sample', 'library'])
             .applymap('{:.1e}'.format)  # scientific notation
             .to_html()
             ))
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>fate</th>
      <th>failed chastity filter</th>
      <th>invalid barcode</th>
      <th>low quality barcode</th>
      <th>unparseable barcode</th>
      <th>valid barcode</th>
    </tr>
    <tr>
      <th>sample</th>
      <th>library</th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th rowspan="2" valign="top">WOO-3-BA.1-FLAG-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>7.4e+05</td>
      <td>2.5e+06</td>
      <td>1.6e+05</td>
      <td>5.7e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>6.5e+05</td>
      <td>2.3e+06</td>
      <td>1.6e+05</td>
      <td>5.3e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-3-BA.1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>8.8e+05</td>
      <td>3.0e+06</td>
      <td>1.9e+05</td>
      <td>6.8e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>7.6e+05</td>
      <td>2.8e+06</td>
      <td>1.9e+05</td>
      <td>6.5e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-3-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>9.9e+05</td>
      <td>1.6e+06</td>
      <td>1.0e+05</td>
      <td>3.3e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>8.6e+05</td>
      <td>1.3e+06</td>
      <td>7.6e+04</td>
      <td>2.6e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-4-BA.1-FLAG-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>5.3e+05</td>
      <td>1.8e+06</td>
      <td>1.1e+05</td>
      <td>4.1e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>4.8e+05</td>
      <td>1.8e+06</td>
      <td>1.2e+05</td>
      <td>4.0e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-4-BA.1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>4.2e+05</td>
      <td>1.5e+06</td>
      <td>9.0e+04</td>
      <td>3.3e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>8.3e+05</td>
      <td>2.8e+06</td>
      <td>1.9e+05</td>
      <td>6.6e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-4-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>9.9e+05</td>
      <td>1.6e+06</td>
      <td>1.1e+05</td>
      <td>3.3e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>1.2e+06</td>
      <td>1.9e+06</td>
      <td>1.1e+05</td>
      <td>3.7e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-5-BA.1-FLAG-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>6.9e+05</td>
      <td>2.4e+06</td>
      <td>1.5e+05</td>
      <td>5.4e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>5.8e+05</td>
      <td>2.1e+06</td>
      <td>1.5e+05</td>
      <td>4.9e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-5-BA.1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>9.3e+05</td>
      <td>3.1e+06</td>
      <td>2.0e+05</td>
      <td>7.2e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>7.0e+05</td>
      <td>2.6e+06</td>
      <td>1.8e+05</td>
      <td>5.9e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-5-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>8.5e+05</td>
      <td>1.4e+06</td>
      <td>9.1e+04</td>
      <td>2.8e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>9.4e+05</td>
      <td>1.4e+06</td>
      <td>8.4e+04</td>
      <td>2.9e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-6-BA.1-FLAG-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>6.3e+05</td>
      <td>2.1e+06</td>
      <td>1.3e+05</td>
      <td>4.9e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>6.6e+05</td>
      <td>2.4e+06</td>
      <td>1.6e+05</td>
      <td>5.5e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-6-BA.1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>8.6e+05</td>
      <td>2.8e+06</td>
      <td>1.8e+05</td>
      <td>6.5e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>9.1e+05</td>
      <td>3.2e+06</td>
      <td>2.2e+05</td>
      <td>7.4e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WOO-6-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>8.7e+05</td>
      <td>1.4e+06</td>
      <td>9.1e+04</td>
      <td>2.8e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>1.1e+06</td>
      <td>1.7e+06</td>
      <td>9.9e+04</td>
      <td>3.3e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WWW-1-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>1.0e+06</td>
      <td>1.7e+06</td>
      <td>1.1e+05</td>
      <td>3.3e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>1.4e+06</td>
      <td>2.1e+06</td>
      <td>1.2e+05</td>
      <td>4.1e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">WWW-2-Wuhan-Hu-1-Strep-abneg</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>8.3e+05</td>
      <td>1.4e+06</td>
      <td>8.8e+04</td>
      <td>2.7e+06</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>1.3e+06</td>
      <td>2.0e+06</td>
      <td>1.1e+05</td>
      <td>4.0e+06</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">none-BA.1-none-ref</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>3.5e+06</td>
      <td>1.1e+07</td>
      <td>7.2e+05</td>
      <td>2.6e+07</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>3.2e+06</td>
      <td>1.1e+07</td>
      <td>7.8e+05</td>
      <td>2.6e+07</td>
    </tr>
    <tr>
      <th rowspan="2" valign="top">none-Wuhan-Hu-1-none-ref</th>
      <th>lib1</th>
      <td>0.0e+00</td>
      <td>7.4e+06</td>
      <td>1.2e+07</td>
      <td>8.2e+05</td>
      <td>2.5e+07</td>
    </tr>
    <tr>
      <th>lib2</th>
      <td>0.0e+00</td>
      <td>1.2e+07</td>
      <td>1.2e+07</td>
      <td>7.2e+05</td>
      <td>1.9e+07</td>
    </tr>
  </tbody>
</table>


Now we plot the barcode-read fates for each library / sample, showing the bars for valid barcodes in orange and the others in gray.
We see that the largest fraction of barcode reads correspond to valid barcodes, and most of the others are invalid barcodes (probably because the map to variants that aren't present in our variant table since we didn't associate all variants with barcodes). The exception to this is lib2 Titeseq_03_bin3; the PCR for this sample in the original sequencing run failed, so we followed it up with a single MiSeq lane. We did not filter out the PhiX reads from this data before parsing, so these PhiX reads will deflate the fraction of valid barcode reads as expected, but does not indicate any problems.


```python
ncol = 4
nfacets = len(fates.groupby(['sample', 'library']))

barcode_fate_plot = (
    ggplot(
        fates
        .assign(sample=lambda x: pd.Categorical(x['sample'],
                                                x['sample'].unique(),
                                                ordered=True),
                fate=lambda x: pd.Categorical(x['fate'],
                                              x['fate'].unique(),
                                              ordered=True),
                is_valid=lambda x: x['fate'] == 'valid barcode'
                ), 
        aes('fate', 'count', fill='is_valid')) +
    geom_bar(stat='identity') +
    facet_wrap('~ sample + library', ncol=ncol) +
    scale_fill_manual(CBPALETTE, guide=False) +
    theme(figure_size=(3.25 * ncol, 2.5 * math.ceil(nfacets / ncol)),
          axis_text_x=element_text(angle=90),
          panel_grid_major_x=element_blank()
          ) +
    scale_y_continuous(labels=dms_variants.utils.latex_sci_not,
                       name='number of reads')
    )

_ = barcode_fate_plot.draw()
```


    
![png](aggregate_variant_counts_files/aggregate_variant_counts_26_0.png)
    


## Initialize codon variant table
For each primary target, do the following:

Initialize the [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) using the wildtype gene sequence and the CSV file with the table of variants:

Now we use the [CodonVariantTable.add_sample_counts_df](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable.add_sample_counts_df) method to add the barcode counts to the variant table:

The variant table now has a `variant_count_df` attribute that gives a data frame of all the variant counts.
Here are the first few lines:

Write the variant counts data frame to a CSV file.
It can then be used to re-initialize a [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) via its [from_variant_count_df](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable.from_variant_count_df) method:


```python
# loop through each target
for target in config['input_codon_variant_tables'].keys():

    wt_seqrecord = Bio.SeqIO.read(config['wildtype_sequence'][target], 'fasta')
    geneseq = str(wt_seqrecord.seq)
    primary_target = wt_seqrecord.name
    print(f"Read sequence of {len(geneseq)} nt for {primary_target} from {config['wildtype_sequence'][target]}")

    print(f"Initializing CodonVariantTable from gene sequence and {config['codon_variant_tables'][target]}")

    variants = dms_variants.codonvarianttable.CodonVariantTable(
                    geneseq=geneseq,
                    barcode_variant_file=config['codon_variant_tables'][target],
                    substitutions_are_codon=True,
                    substitutions_col='codon_substitutions',
                    primary_target=primary_target)

    variants.add_sample_counts_df(counts)

    display(HTML(variants.variant_count_df.head().to_html(index=False)))

    print(f"Writing variant counts to {config['variant_counts'][target]}")
    variants.variant_count_df.to_csv(config['variant_counts'][target], index=False, compression='gzip')
```

    Read sequence of 603 nt for Wuhan-Hu-1 from data/Wuhan-Hu-1_sequence.fasta
    Initializing CodonVariantTable from gene sequence and results/variants/codon_variant_table_Wuhan_Hu_1.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>sample</th>
      <th>barcode</th>
      <th>count</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>GAAACAAATTTCTATA</td>
      <td>5871</td>
      <td>16</td>
      <td>GTA52GGT</td>
      <td>V52G</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>GACTATCGAATTATTG</td>
      <td>5120</td>
      <td>10</td>
      <td>AAG48TCT</td>
      <td>K48S</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>TTAAAAGCAACCGACC</td>
      <td>5078</td>
      <td>12</td>
      <td>TCT36CAT</td>
      <td>S36H</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>CGTAACATTTACATAT</td>
      <td>4863</td>
      <td>14</td>
      <td>ACT15TAT</td>
      <td>T15Y</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>AATGAACCCAAAACAA</td>
      <td>4775</td>
      <td>18</td>
      <td>AGA16TTT</td>
      <td>R16F</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


    Writing variant counts to results/counts/variant_counts_Wuhan_Hu_1.csv.gz
    Read sequence of 603 nt for BA.1 from data/BA1_sequence.fasta
    Initializing CodonVariantTable from gene sequence and results/variants/codon_variant_table_BA1.csv



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>sample</th>
      <th>barcode</th>
      <th>count</th>
      <th>variant_call_support</th>
      <th>codon_substitutions</th>
      <th>aa_substitutions</th>
      <th>n_codon_substitutions</th>
      <th>n_aa_substitutions</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>BA.1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>ACAATCCATGCACTAA</td>
      <td>1233</td>
      <td>6</td>
      <td>TTG38AGA</td>
      <td>L38R</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA.1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>TATACTAAACACATAA</td>
      <td>423</td>
      <td>3</td>
      <td>TTC160AGA</td>
      <td>F160R</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA.1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>AATTAAAAATAAAAAT</td>
      <td>332</td>
      <td>23</td>
      <td>ACT201TGT</td>
      <td>T201C</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA.1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>CATGCAAAGACCTCGA</td>
      <td>197</td>
      <td>4</td>
      <td>ATA142AAT</td>
      <td>I142N</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>BA.1</td>
      <td>lib1</td>
      <td>none-Wuhan-Hu-1-none-ref</td>
      <td>CAACCCACAAACTCAA</td>
      <td>193</td>
      <td>1</td>
      <td>GAT34ACT</td>
      <td>D34T</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


    Writing variant counts to results/counts/variant_counts_BA1.csv.gz


The [CodonVariantTable](https://jbloomlab.github.io/dms_variants/dms_variants.codonvarianttable.html#dms_variants.codonvarianttable.CodonVariantTable) has lots of nice functions that can be used to analyze the counts it contains.
However, we do that in the next notebook so we don't have to re-run this entire (rather computationally intensive) notebook every time we want to analyze a new aspect of the counts.
