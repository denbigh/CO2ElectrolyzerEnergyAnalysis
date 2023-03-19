This folder contains code for generation of all plots in "Electrolyzer Energy 
Dominates Separation Costs in State-of-the-Art CO2 Electrolyzers: Implications 
for Single-Pass CO2 Utilization" by Moore et al. 

# Electrolyzer_Model.py
Process model. Implements electrolyzer and gas separation calculations, global 
mass and energy balances, and TEA calculations. Run using Python kernel 3.7.11.  

# PaperFigs.ipynb
A Jupyter notebook which generates the majority of plots in the manuscript and 
SI. Most Figures utilize models in Electrolyzer_Model.py. Run using Python 
kernel 3.7.11.  

# KimetalImplementation.ipynb
A Jupyter notebook which calculates stripper heat duty via the methodology of 
Kim et al., and generates Figures 5, S6, S7 and S8. Run using Julia kernel 
1.7.1.

# License
Copyright (c) 2023 Thomas Moore

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
