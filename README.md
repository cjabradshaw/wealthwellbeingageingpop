# Effects of ageing/declining populations on national-level indices of socio-economic performance
<a href="https://doi.org/10.5281/zenodo.15826278"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.15826278.svg" alt="DOI"></a>
<img align="right" src="www/healthwellbeing.png" alt="health, wealth, wellbeing" width="250" style="margin-top: 20px">

Corey J. A. Bradshaw<br>
<a href="http://globalecologyflinders.com">Global Ecology</a>, Flinders University, Australia<br>
<a href="mailto:corey.bradshaw@flinders.edu.au">e-mail</a><br>
<br>
Shana M. McDermott<br>
<a href="https://economics.missouri.edu">Department of Economics</a>, University of Missouri, USA<br>
<a href="mailto:smmvt@missouri.edu">e-mail</a><br>
<br>
Matthew E. Oliver<br>
<a href="https://econ.gatech.edu">School of Economics</a>, <a href="https://www.gatech.edu">Georgia Institute of Technology</a>, USA<br>
<a href="mailto:matthew.oliver@econ.gatech.edu">e-mail</a><br>
<br>
Accompanies pre-print:<br>
Bradshaw, CJA, SM McDermott. <a href="https://doi.org/10.48550/arXiv.2508.16872">No evidence ageing or declining populations compromise socio-economic performance of countries</a>. <em>arχiv</em> doi:10.48550/arXiv.2508.16872<br>

## <a href="https://github.com/cjabradshaw/wealthwellbeingageingpop/tree/main/scripts">Scripts</a>
- <code>ageingSEperformanceGH.R</code>: R code for all analyses
 
## <a href="https://github.com/cjabradshaw/wealthwellbeingageingpop/tree/main/data">Data</a>
- <em>continent.country2.csv</em>: country names, 3-character ISO country codes, continental region, other regionalisation
- <em>CPIA23wb.csv</em>: CPIA transparency, accountability, and corruption in the public sector rating (source: World Bank)
- <em>cpi.csv</em>: corruption perception index (source: Transparency International)
- <em>cpits.csv</em>: corruption perception index time series (source: Transparency International)
- <em>DCWI.csv</em>: per-capita domestic comprehensive wealth index by country (source: World Bank)
- <em>freedom.csv</em>: freedom score time series (source: Freedom House)
- <em>freedom2025.csv</em>: freedom score 2025 (source: Freedom House)
- <em>gdppcPPP.csv</em>: per-capita gross domestic product adjusted for purchasing power parity by country (source: World Bank)
- <em>giniMn.csv</em>: annual time series of Gini coefficient (income equality) (source: World Bank)
- <em>gcf.csv</em>: gross captial formation (source: World Bank)
- <em>gor.csv</em>: grants and other revenue (source: World Bank)
- <em>HDI.csv</em>: Human Development Index (source: United Nations Development Programme)
- <em>HDIPP.csv</em>: planetary pressure-adjusted Human Development Index (source: United Nations Development Programme)
- <em>HDIPPtimeseries.csv</em>: planetary pressure-adjusted Human Development Index time series (source: United Nations Development Programme)
- <em>healthyLE.csv</em>: healthy life expectancy at birth (years) (source: World Health Organization)
- <em>mva.csv</em>: value-added manufacturing (source: World Bank)
- <em>par.csv</em>: patent applications (source: World Bank)
- <em>polspectrum.csv</em>: political spectrum (Global Parliament Index 2025) (source: Arden Strategies)
- <em>pwt.csv</em>: human capital index, capital services, capital stock, output-side real GDP at current PPP, capital services levels at current PPPs (source: Penn World Table, Groningen Growth and Development Centre, Faculty of Economics and Business)
- <em>popXage.csv</em>: population size by country, year (1950-2021), and yearly age class (0-100+) (source: United Nations Population Division)
- <em>rde.csv</em>: gross expenditure on research and development (source: World Bank)
- <em>savings.csv</em>: gross savings (source: World Bank)
- <em>WB_ASPD.csv.zip</em>: total factor productivity (zipped) (source: World Bank)
- <em>WB_ASPD_LPXR.csv</em>: labour productivity (source: World Bank)
- <em>WB_WDI-SI_POV_GINI.csv</em>: Gini coefficient (income equality) (source: World Bank)
- <em>wellbeingrank.csv</em>: composite wellbeing rank by country (source: Blanchflower & Bryson 2024)
 
## Required R libraries
<code>boot</code>, <code>car</code>, <code>dismo</code>, <code>gbm</code>, <code>ggplot2</code>, <code>ggpubr</code>, <code>ggrepel</code>, <code>ggtext</code>, <code>lmtest</code>, <code>nlme</code>, <code>sandwich</code>, <code>usdm</code>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University logo" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="130" style="margin-top: 20px"></a>  &nbsp; &nbsp;
 <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHF_Logo_Email_Version Transparent.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; &nbsp; <a href="https://economics.missouri.edu"><img align="bottom-left" src="www/UMlogo.png" alt="U Missouri logo" width="180" style="margin-top: 20px"></a></p>
