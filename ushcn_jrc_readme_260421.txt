Data files to accompany the analysis of summer and winter extremes published here: https://link.springer.com/article/10.1007/s00704-026-06200-3

There are three data files in this directory.

ushcn_station_list_250617.txt

This is the list of 1,218 United States Historical Climate Network stations by numerical order of their heritage Cooperative Observer Program (COOP) station identifier.  This identifier has 6 numbers, the first two being the state identifier and the last four being the individual station identifier, usually corresponding to alphabetical order of the station name.  However, note that the first nine states (AL to GA) display only five COOP ID digits as the leading “0” does not appear in numerically-based files such as this.  Thus “11084” in this file is “011084” in the official USHCN station listing. 

ushcn_jrc_tmin_260421.txt	[Dec, Jan, Feb, Mar of each year b. Dec 1898]
ushcn_jrc_tmax_260421.txt	[May, Jun, Jul, Aug, Sep of each year b. May 1899]	

In these files are the values of daily minimum (tmin) or maximum (tmax) temperature for each station with one month per line in integer degrees Fahrenheit which is the native observed metric.  “-999” is “missing”.   The filename includes the date of the latest version.  The filename includes the date the data were provided (260421 is 21 Apr 2026).  Revisions are possible as new data may be added to fill-in some of the remaining gaps.  The data are updated through 30 Sept 2025.  The next update will include data through 30 Sep 2026 which won't be available until late 2026 or early 2027.   

Each record consists of NNNNNN  YYYY  MM [followed by 31 values] NTEMP
Where:
NNNNNN = Station COOP ID
YYYY = Year
MM = Month
NTEMP(31) = 31 values of temperature for the specific month.

Note that the two datafiles are not in numerical order by COOP station ID as are found in the ushcn_jrc_station_list.txt file but are in the order that NOAA/NCEI provides when downloading (see below) in which those stations for which NCEI switched to WBAN identifiers (about 80) are at the end. An example of both the WBAN change and missing leading “0” is for Fresno Yosemite AP CA which is listed by its COOP ID in all files here as 43257 (“04” is California), but occurs near the end of the temperature datafiles because the NCEI station order uses the WBAN ID of 93193.

Values are downloaded from NOAA/NCEI for the original 1,218 stations here 
ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/daily/hcn/ and described in general here 
https://www.ncei.noaa.gov/products/land-based-station/us-historical-climatology-network

Please note we use the daily datafiles, not the monthly homogenized files (52j version).

NCEI has also identified over 200 stations that are compatible and are used to fill gaps or extend the record back in time for appropriation stations.  These additional data are available from sources such as XMACiS and were merged or “threaded” with the original USHCN data records.  These two groups of datafiles account for about 90 percent of the values in ushcn_jrc datafiles.  In rare cases several months of the original USHCN data were set to missing and in-filled with nearby bias-corrected station values due to suddenly sporadic reporting or poor values.

Where there were still gaps, nearby stations with high correlation were accessed and merged with the existing data – in essence following NOAA/NCEI threading procedure.  Several tens of thousands of these station values were manually keyed-in as they had not been digitized.  The database consists of daily minimum temperatures in the cold season (Dec to Mar) and daily maximum temperatures in the warm season (May to Sep).  

One main purpose of this dataset is to discover information about extreme values rather than to calculate small trends.  It is well known that the minimum values have been affected (warmed) considerably by the growing presence of manufactured surfaces around the stations (concrete, asphalt, buildings, etc.) and this is borne out in the results.  Maximum temperatures tend to be less influenced by these changes as the daytime mixing of air both vertically and horizontally minimizes this human impact.  There are also discontinuities within stations due to relocations, or changes in instruments or time of observation.  These changes are considered small (up to 1 to 2 °F) compared to the magnitude of the extremes we seek to identify (10’s of °F).

The most common situation requiring supplemental data occurred with the closure of a USHCN station.  Values of a nearby, highly correlated station were merged with the closed station once a bias was determined and removed (separately for Tmin and Tmax).  The supplemental station was often an airport station or in the case of Oklahoma, an Oklahoma Mesonet station.  

In many of the analyses which utilized this database, a requirement of 92 percent of the data was set, providing sufficient data records for 1,211 of the 1,218 stations.  Three stations in West TX, two in NY (data yet to be manually keyed) and one each in WA and UT do not have 92 percent of the observations.  Unlike NY, the other five stations are likely not producible because of lack of credible nearby observations back to 1898.  

Death Valley (Greenland Ranch) data for several years around 1913 were subsitituted from nearby stations due to obvious errors (see https://journals.ametsoc.org/view/journals/bams/aop/BAMS-D-24-0313.1/BAMS-D-24-0313.1.pdf)

In constructing the database, it was required that (a) data for 1898-1899 be present to capture the Feb 1899 Arctic Outbreak, and (b) data for 2021-2025 be present to place recent extremes in proper perspective.  The median value of available data per station is 98 percent.

With on the order 40,000,000 temperature values there are certainly erroneous entries.  We have discovered and corrected many, such as true hot temperature spikes in the summer along the Pacific Coast that were erroneously set to missing by the NCEI algorithms.  Another example is that of many stations in the Southern Plains during Feb 1899 when those keying in the data misread the notations of the observers who intended to provide below zero values, but invented their own systems, such as “05” which to the observer meant “5 below zero” but was keyed in as +5.    

It would be greatly appreciated that if questionable values are discovered that these could be brought to our attention so that we may correct or eliminate them.  Send information to christy@nsstc.uah.edu with subject line “USHCN.”  Also, any questions about the dataset may be sent to the same address.

An example of the FORTRAN code used to calculate the heat/cold waves and some of the percentile statistics is in us_dly_waves.f.  Because no adjustments were applied, some of the code is superfluous. It is virtually certain the user will have his/her own code (much more efficient than FORTRAN) to process the datafiles as needed so that machine-specific data will in this code will be unnecessary.  

John R Christy
The University of Alabama in Huntsville
21 April 2026


