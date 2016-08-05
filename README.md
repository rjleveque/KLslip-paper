# KLslip-paper

Python code to accompany the paper

- Generating Random Earthquake Events for Probabilistic Tsunami Hazard Assessment
  by R. J. LeVeque, K. Waagan, F. I. Gonz ́alez, D. Rim, and G. Lin,
  To appear in Pure & Applied Geophysics.

The scripts `KL1d_figures.py` and `KL2d_figures.py` produce the figures that
appear in the paper. 

The Jupyter notebooks `KL1d_figures_demo.ipynb` and `KL2d_figures_demo.ipynb`
contain essentially the same code in a more readable form 
with some additional comments.  Also the notebook
`KL2d_figures_truncations.ipynb` was used to produce Figure 13.

These notebooks were committed to the repository with
figures included so you can view them on the GitHub pages:

- https://github.com/rjleveque/KLslip-paper/blob/master/plot_chile2010_demo.ipynb
- https://github.com/rjleveque/KLslip-paper/blob/master/KL1d_figures_demo.ipynb
- https://github.com/rjleveque/KLslip-paper/blob/master/KL2d_figures_demo.ipynb
- https://github.com/rjleveque/KLslip-paper/blob/master/KL2d_figures_truncations.ipynb 



To run these, the following dependencies are needed:

- Python 2.x with the packages:
  Numpy, Scipy, matplotlib, pandas

- Clawpack Version 5.3.1, see http://www.clawpack.org

- The file `dtopotools.py` in this repository is a version modified from
  what is in GeoClaw in Version 5.3.1.  
  The figures in the paper were generated using commit e21399557b8c537.

Note: The seed is set for the random number generator in order to reproduce
the figures exactly.  Changing the seed will produce different realizations
and may change the plots.

