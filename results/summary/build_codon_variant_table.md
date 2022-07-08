# Build combined 'codon_variant_table.csv'

Need to take the codon variant tables from the "epistatic shifts" _Science_ paper and from the unpublished Omicron (BA.1 and BA.2) projects, and extract and combine the barcodes that correspond to these variants:
* Wuhan-Hu-1
* BA.1

Also need to duplicate the "follow-up pools" for Wuhan-Hu-1 (if I can find where those are) from the _Science_ paper and duplicate into library 2. 

The problem is this: it is not clear which barcodes came from the follow-up pool, and which did not. 
The follow-up pool was sequenced separately, so might be able to use this to help.

Ideal would be to find barcodes that correspond to the follow-up pool, query the codon variant table on those barcodes, duplicate them but assign to library 2. 

I'm not going to worry about this right now, and will just focus on what we already have.


```python
import os

import numpy
import pandas as pd
from IPython.display import display, HTML

import yaml
```

Read config file.


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Define input and output directories.


```python
datadir = 'data'
resultsdir = config['variants_dir']

os.makedirs(resultsdir, exist_ok=True)
os.makedirs(config['prior_DMS_data_dir'], exist_ok=True)
```

Read in the new input codon variant tables and assign a target for each.


```python
# hard-coded dictionary to replace values in codon variant table:

names_df = {'Wuhan_Hu_1':'Wuhan-Hu-1',
            'BA1':'BA.1',
            'pool1':'lib1',
            'pool2':'lib2',
           }
```


```python
codon_variant_table = pd.DataFrame()

for target, file in config['input_codon_variant_tables'].items():
    print(f'Reading codon variant table for {target}.')
    df = pd.read_csv(file).replace(names_df)
    codon_variant_table = pd.concat([codon_variant_table, df])
    df.to_csv(config['codon_variant_tables'][target], index=False)

codon_variant_table.to_csv(config['codon_variant_table'], index=False)
display(HTML(codon_variant_table.head().to_html(index=False)))
```

    Reading codon variant table for Wuhan-Hu-1.
    Reading codon variant table for BA.1.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>library</th>
      <th>barcode</th>
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
      <td>AAAAAAAAAAAGGAGA</td>
      <td>4</td>
      <td>GGT166ATG</td>
      <td>G166M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>AAAAAAAAAAATTTAA</td>
      <td>4</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>0</td>
      <td>0</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>AAAAAAAAAACGCGTA</td>
      <td>3</td>
      <td>GAA154ACT</td>
      <td>E154T</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>AAAAAAAAAACTCCAA</td>
      <td>2</td>
      <td>TTT156ATG</td>
      <td>F156M</td>
      <td>1</td>
      <td>1</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>lib1</td>
      <td>AAAAAAAAACCGATTA</td>
      <td>2</td>
      <td>CAG84GAA</td>
      <td>Q84E</td>
      <td>1</td>
      <td>1</td>
    </tr>
  </tbody>
</table>


## Now do the same for the DMS ACE2 binding and RBD expression scores for each library. 


```python
# hard-coded dictionary to replace values in final variant scores:

names_df = {'Wuhan_Hu_1':'Wuhan-Hu-1',
            'Omicron_BA1':'BA.1',
           }
```


```python
mut_bind_expr = pd.DataFrame()

for target, file in config['input_mut_bind_expr'].items():
    print(f'Reading variant scores for {target}.')
    df = (pd.read_csv(file)
          .replace(names_df)
          .query('target==@target')
         )
    mut_bind_expr = pd.concat([mut_bind_expr, df])

mut_bind_expr.to_csv(config['mut_bind_expr'], index=False)
display(HTML(mut_bind_expr.head().to_html(index=False)))
```

    Reading variant scores for Wuhan-Hu-1.
    Reading variant scores for BA.1.



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>target</th>
      <th>wildtype</th>
      <th>position</th>
      <th>mutant</th>
      <th>mutation</th>
      <th>bind</th>
      <th>delta_bind</th>
      <th>n_bc_bind</th>
      <th>n_libs_bind</th>
      <th>bind_rep1</th>
      <th>bind_rep2</th>
      <th>bind_rep3</th>
      <th>expr</th>
      <th>delta_expr</th>
      <th>n_bc_expr</th>
      <th>n_libs_expr</th>
      <th>expr_rep1</th>
      <th>expr_rep2</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>N</td>
      <td>331</td>
      <td>A</td>
      <td>N331A</td>
      <td>8.79360</td>
      <td>0.06027</td>
      <td>4</td>
      <td>2</td>
      <td>8.76603</td>
      <td>NaN</td>
      <td>8.82117</td>
      <td>10.29895</td>
      <td>0.11422</td>
      <td>2</td>
      <td>1</td>
      <td>10.29895</td>
      <td>NaN</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>N</td>
      <td>331</td>
      <td>C</td>
      <td>N331C</td>
      <td>8.61594</td>
      <td>-0.15567</td>
      <td>5</td>
      <td>3</td>
      <td>8.73710</td>
      <td>8.56255</td>
      <td>8.54816</td>
      <td>9.67665</td>
      <td>-0.50923</td>
      <td>4</td>
      <td>2</td>
      <td>9.49750</td>
      <td>9.85579</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>N</td>
      <td>331</td>
      <td>D</td>
      <td>N331D</td>
      <td>8.75409</td>
      <td>-0.01751</td>
      <td>8</td>
      <td>3</td>
      <td>8.65990</td>
      <td>8.79668</td>
      <td>8.80570</td>
      <td>10.06985</td>
      <td>-0.11602</td>
      <td>5</td>
      <td>2</td>
      <td>10.14610</td>
      <td>9.99361</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>N</td>
      <td>331</td>
      <td>E</td>
      <td>N331E</td>
      <td>8.92561</td>
      <td>0.15400</td>
      <td>10</td>
      <td>3</td>
      <td>8.69116</td>
      <td>9.12888</td>
      <td>8.95680</td>
      <td>10.18436</td>
      <td>-0.00151</td>
      <td>6</td>
      <td>2</td>
      <td>10.22575</td>
      <td>10.14298</td>
    </tr>
    <tr>
      <td>Wuhan-Hu-1</td>
      <td>N</td>
      <td>331</td>
      <td>F</td>
      <td>N331F</td>
      <td>8.65690</td>
      <td>-0.11470</td>
      <td>6</td>
      <td>3</td>
      <td>8.36984</td>
      <td>8.80036</td>
      <td>8.80051</td>
      <td>10.01397</td>
      <td>-0.17191</td>
      <td>4</td>
      <td>2</td>
      <td>10.14360</td>
      <td>9.88434</td>
    </tr>
  </tbody>
</table>



```python

```
