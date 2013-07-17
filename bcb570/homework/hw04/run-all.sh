#!/bin/bash
curl -o wdbc.data http://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/wdbc.data
perl -ne '@f = split/,/; print(join(",", @f[2..31]))' < wdbc.data > wdbc-values.data
perl -ne '@f = split/,/; printf("%s\n", $f[1])' < wdbc.data > wdbc-diagnoses.data
/Applications/MATLAB_R2011b.app/bin/matlab -nodisplay < run-kmeans.m
echo -e "\n\n\n======Results======"
paste -d: wdbc-clusters.data wdbc-diagnoses.data | sort | uniq -c
