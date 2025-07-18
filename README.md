# Effects of ageing/declining populations on national-level indices of wealth and wellbeing
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

## <a href="https://github.com/cjabradshaw/wealthwellbeingageingpop/tree/main/scripts">Scripts</a>
- <code>wealthwellbeingpopGH.R</code>: R code for all analyses
 
## <a href="https://github.com/cjabradshaw/wealthwellbeingageingpop/tree/main/data">Data</a>
- <em>continent.country2.csv</em>: country names, 3-character ISO country codes, continental region, other regionalisation
- <em>DCWI.csv</em>: per-capita domestic comprehensive wealth index by country (source: World Bank)
- <em>gdppcPPP.csv</em>: per-capita gross domestic product adjusted for purchasing power parity by country (source: World Bank)
- <em>HDI.csv</em>: Human Development Index (source: United Nations Development Programme)
- <em>HDIPP.csv</em>: planetary pressure-adjusted Human Development Index (source: United Nations Development Programme)
- <em>popXage.csv</em>: population size by country, year (1950-2021), and yearly age class (0-100+) (source: United Nations Population Division)
- <em>wellbeingrank.csv</em>: composite wellbeing rank by country (source: Blanchflower & Bryson 2024)
 
## Required R libraries
<code>boot</code>, <code>dismo</code>, <code>gbm</code>, <code>ggarrange</code>, <code>ggplot2</code>, <code>ggpubr</code>, <code>ggrepel</code>, <code>usdm</code>

<p><a href="https://www.flinders.edu.au"><img align="bottom-left" src="www/Flinders_University_Logo_Stacked_RGB_Master.jpg" alt="Flinders University logo" width="80" style="margin-top: 20px"></a> &nbsp; <a href="https://globalecologyflinders.com"><img align="bottom-left" src="www/GEL Logo Kaurna New Transp.png" alt="GEL logo" width="130" style="margin-top: 20px"></a>  &nbsp; &nbsp;
 <a href="https://ciehf.au"><img align="bottom-left" src="www/CIEHF_Logo_Email_Version Transparent.png" alt="CIEHF logo" width="200" style="margin-top: 20px"></a>  &nbsp; &nbsp; &nbsp; <a href="https://economics.missouri.edu"><img align="bottom-left" src="www/UMlogo.png" alt="U Missouri logo" width="180" style="margin-top: 20px"></a></p>
