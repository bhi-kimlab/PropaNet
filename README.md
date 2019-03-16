# PropaNet
This program was designed for...

## Installation
To download all the examples, simply clone this repository
`https://github.com/minwoopak/Propanet.git`

## Dependency
To run them, you also need the following dependencies

---
#### Python
 * numpy
 * pandas
 * multiprocessing

#### R
 * Limma

## Input File Format
Time-series gene expression data of multiple conditions has to be stored in a single tab-delimited matrix form. (where different conditions are stacked column-wise)

| Gene_id | Control          | Condition        |
| ------- | :--------------: | ---------------: |
| Gene1   | expression level | expression level |
| Gene2   | expression level | expression level |
| ...     | ...              | ...              |

### Sample Labels
The input file contains a column label as its first line. The label has a format of i_j_h
 * i : Control/Condition
 * j : time-point
 * h : time-unit (hour)

### Gene IDs
The first column of the file includes gene IDs.

### Example
```
AGI Control_0_h Control_05_h    Control_1_h Control_3_h Control_6_h Control_12_h    Control_24_h    Cold_0_h    Cold_05_h   Cold_1_h    Cold_3_h    Cold_6_h    Cold_12_h   Cold_24_h
AT1G01010   0   -0.139394407    -0.152530099    -0.156121134    -0.180087384    -0.095033671    -0.216560433    0   -0.154829544    -0.070773114    -0.146864574    -0.142279603    -0.160496473    -0.110211809
AT1G01040   0   -0.042997223    0.059865445 -0.030210557    -0.07748179 -0.045007578    0.021080998 0   -0.045927327    -0.020536712    0.015591588 0.064597276 0.067403396 0.136353683
AT1G01060   0   -0.05952888 -0.083900109    -0.417545241    -0.723319018    -0.691957548    -0.01210007 0   -0.042703729    -0.000666969    -0.02564026 -0.034619052    -0.074636705    -0.075214348
AT1G01070   0   -0.016709256    -0.071140518    -0.010718778    -0.295093989    -0.283263794    -0.049080058    0   0.016980682 -0.036470699    -0.05335869 -0.137604399    -0.212629356    -0.245722342
AT1G01080   0   -0.008805295    0.005914988 0.014547664 0.027181841 0.081628037 0.020220747 0   0.035903529 -0.009913287    -0.0011389  -0.009907133    -0.049022001    -0.084337305
```

---
## Run Propanet
```bash
bash Step0_DEG_extraction.sh
python TF_adding_NP.py
bash run_TF.sh
```
