#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Taken from https://github.com/supl/xlsx2tsv
# Modified: 2018-01-25 10:57:28 - data_only to get values and not formulas

import csv
import xlrd
import sys

if len(sys.argv) < 2: sys.exit(-1)

filename = sys.argv[1]
workbook = xlrd.load_workbook(filename,data_only=True) #data_only to get values and not formulas
sheetname = workbook.get_sheet_names()[0] # retrive first sheet by defaul

worksheet = workbook.get_sheet_by_name(sheetname)
csvwriter = csv.writer(sys.stdout, delimiter="\t")
for row in worksheet.rows:
    stringized = map(lambda x: str(x.value) if x.value is not None else "", row)
    csvwriter.writerow(stringized)
