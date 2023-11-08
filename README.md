# Tropical-eastern-Atlantic-Benthic-Foraminifera-analyses
These scripts were used to process the benthic foraminifera (BF) data GeoB9512-5 of the tropical eastern Atlantic (TE Atlantic).
This folder contains 2 codes, as a supplement to the publication xxxx:

# Folder 9512Diversity_NMDS

_1. Code 9512_MultivariateBFAnalisys.R_
To run and visualize the NMDS of benthic foraminifera from the last deglaciation identified onsite GeoB9512-5

- 9512Sp data comes from Barragán-Montilla, S. (dataset in review b). Benthic Foraminifera counts off NW Africa during the last deglaciation. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.962951
- 9512Gr2.csv as reported in the Supplementary Information
- 9512EA.csv adapted for this code using the Supplementary Information 3 y 4

# Folder 9512Age_Proxies

_2. Code 9512AgeModel.R was written in 3 parts._

**First Part:** lines 9 - 290 to model ages based on the radiocarbon data published here https://doi.org/10.1594/PANGAEA.962899 
and estimate median ages of the benthic foraminifera geochemistry data points reported here https://doi.org/10.1594/PANGAEA.962968 
with the corresponding uncertainties calculated for this study

**Second part:** from lines 291 - 965, to visualize the results presented in the publication
9512_EA.csv file adapted for this code using Supplementary Information 3 y 4

**A third optional** (lines 968 - 1,115 ) part that includes global variables compiled by different authors, including:

Atlantic Ph/Th Line 964 https://www.nature.com/articles/nature02494 and https://www.nature.com/articles/nature14059
Global CO2 https://doi.pangaea.de/10.1594/PANGAEA.871265 line 1047
Global Mean Seasurface Temprature (GMST) https://www.nature.com/articles/s41586-021-03984-4 line 1078

All necessary files are provided in the folders
Any inquiries can be directed to the corresponding author

## References

- Barragán-Montilla, S. (dataset in review b). Benthic Foraminifera counts off NW Africa during the last deglaciation. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.962951
- Barragán-Montilla, S. (dataset in review c). Benthic foraminifera stable isotopes (d18O and d13C) of sediment core GeoB9512-5. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.962968
- Barragán-Montilla, S. and Mulitza, S. (dataset in review). Mulitza, Stefan: Radiocarbon ages of sediment core GeoB9512-1. PANGAEA, https://doi.pangaea.de/10.1594/PANGAEA.962899
- Blaauw, M. & Christen, J. A. Flexible paleoclimate age-depth models using an autoregressive gamma process. Bayesian Anal. 6, 457–474 (2011).
- Böhm, E. et al. Strong and deep Atlantic meridional overturning circulation during the last glacial cycle. Nature 517, 73–76 (2015).
- Heaton, T. J. et al. Marine20—The Marine Radiocarbon Age Calibration Curve (0–55,000 cal BP). Radiocarbon 62, 779–820 (2020).
- Kranner, M., Harzhauser, M., Beer, C., Auer, G. & Piller, W. E. Calculating dissolved marine oxygen values based on an enhanced Benthic Foraminifera Oxygen Index. Sci Rep 12, 1376 (2022).
- Köhler, P., Nehrbass-Ahles, C., Schmitt, J., Stocker, T.F., & Fischer, H. (2017). A 156kyr smoothed history of the atmospheric greenhouse gases CO2, CH4, and N2O and their radiative forcing. Earth System Science Data, 9, 363–387, https://doi.org/10.5194/essd-9-363-2017.
- McManus, J. F., Francois, R., Gherardi, J.-M., Keigwin, L. D. & Brown-Leger, S. Collapse and rapid resumption of Atlantic meridional circulation linked to deglacial climate changes. Nature 428, 834–837 (2004).
- Mulitza, Stefan; Bickert, Torsten; Bostock, Helen C; Chiessi, Cristiano Mazur; Donner, Barbara; Govin, Aline; Harada, Naomi; Huang, Enqing; Johnstone, Heather J H; Kuhnert, Henning; Langner, Michael; Lamy, Frank; Lembke-Jene, Lester; Lisiecki, Lorraine E; Lynch-Stieglitz, Jean; Max, Lars; Mohtadi, Mahyar; Mollenhauer, Gesine; Muglia, Juan; Nürnberg, Dirk; Paul, André; Rühlemann, Carsten; Repschläger, Janne; Saraswat, Rajeev; Schmittner, Andreas; Sikes, Elisabeth L; Spielhagen, Robert F; Tiedemann, Ralf (2021): World Atlas of late Quaternary Foraminiferal Oxygen and Carbon Isotope Ratios (WA_Foraminiferal_Isotopes_2022). PANGAEA, https://doi.org/10.1594/PANGAEA.936747 
- Osman, M.B., Tierney, J.E., Zhu, J., Tardif, R., Hakim, G. J., King, J., & Poulsen, C.J. (2021). Globally resolved surface temperatures since the Last Glacial Maximum. Nature 599, 239–244. https://doi.org/10.1038/s41586-021-03984-4
- Völpel, R., Mulitza, S., Paul, A., Lynch‐Stieglitz, J. & Schulz, M. Water Mass Versus Sea Level Effects on Benthic Foraminiferal Oxygen Isotope Ratios in the Atlantic Ocean During the LGM. Paleoceanog and Paleoclimatol 34, 98–121 (2019).
- Waelbroeck, C. et al. The timing of the last deglaciation in North Atlantic climate records. Nature 412, 724–727 (2001).

