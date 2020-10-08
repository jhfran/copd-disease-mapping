# Disease mapping using standardised mortality/morbidity ratios (SMR)

This was an MSc assignment with the purpose of investigating England's hospital admission rates for COPD (2001-2010).

## Objectives: 

- Summarise number of hospital admissions for COPD in England between 2001 and 2010
- Estimate raw and smoothed standardised morbidity ratios (SMR) and look for any spatial patterns 
- Investigate any changes in the risks of hospitalisation due to COPD in England over time
- Consider this scenorio:a government minister wants to increase resources to hospitals to cope with admissions for COPD, recommend key areas/local authorities in England where these resources should be allocated  

## Files 

- R file with my analysis
- copdobserved.csv and copdexpected.csv are the datasets used for this project
    * copdobserved.csv has the variables name - Name of local authority, 20XX - Observed number of hospital admissions for COPD in the year 20XX
    * copdexpected.csv has the variables name - Name of local authority, E20XX - Expected number of hospital admissions for COPD in the year 20XX
    * expected number of cases, previously calculated using indirect standardisation by applying the age–sex dpecific rates for the whole of England to the age–sex population profile of each of the local authorities
- shapefile folder containing shapefiles for England split by local authorities provided in order to produce any maps for this analysis

### Report containing recommendations (within Conclusion/Recommendations section) can be found on my github page : https://jhfran.github.io/copd
