
======================================================================
The 2011 Tohoku Earthquake Tsunami Joint Survey Group
======================================================================

At 14:46 local time on March 11, 2011, a magnitude 9.0 earthquake occurred off the coast of northeast Japan. This earthquake caused a tsunami which attacked Japan as well as a wide range of localities around the Pacific Ocean.  Tsunami surveys were conducted by joint research groups with the participation of 299 tsunami, coastal, seismology, and geology researchers from 64 universities and institutes throughout Japan 

The 2011 Tohoku Earthquake Tsunami Joint Survey Group (hereinafter, the survey group) is an autonomous survey organization and it consists of members from different fields of natural science, tsunami engineering, coastal engineering, and tsunami-related research.  The survey group was managed with researchers at the Faculty of Safety Science of Kansai University and the Disaster Prevention Research Institute of Kyoto University (denoted survey secretariat).

Host web site: http://www.coastal.jp/tsunami2011/

======================================================================
 About this dataset
======================================================================

The post-event tsunami survey data were collected by scientific volunteers following the Tohoku earthquake tsunami. The individual data were measured under the severe conditions after the disaster, so please adhere to the following terms, conditions, and general guidelines:

(1) The data are open for public use on the condition that users site   
	a)  "The 2011 Tohoku Earthquake Tsunami Joint Survey Group",  
	b)  web site address, and 
	c)  the release date. 
at least. For example,
"Data are from the 2011 Tohoku Earthquake Tsunami Joint Survey Group, release 20120330, http://www.coastal.jp/ttjt/"

The release date is important because we continuously check and update the data set since May 2011.

(2) Although we paid careful attention to the quality of the data, the quality can vary with location and surveyor. Please check degree of quality of each data.  Inconsistencies should be reported to the contact provided below.

References
- A prompt report from survey group
-- The 2011 Tohoku Earthquake Tsunami Joint Survey Group (2011) Nationwide Field Survey of the 2011 Off the Pacific Coast of Tohoku Earthquake Tsunami, Journal of Japan Society of Civil Engineers, Series B, Vol. 67 (2011) , No.1 pp.63-66.
- Short scientific summary of survey
-- Mori, N., T.Takahashi, T.Yasuda and H.Yanagisawa (2011) Survey of 2011 Tohoku earthquake tsunami inundation and run-up, Geophysical Research Letters, 38, L00G14, doi:10.1029/2011GL049210.
--Detail analysis of survey dataset
-- Nobuhito Mori, Tomoyuki Takahashi and The 2011 Tohoku Earthquake Tsunami Joint Survey Group (2012) Nationwide survey of the 2011 Tohoku earthquake tsunami, Coastal Engineering Journal, in press, Vol.54, Issue 1, pp.1-27.  (doi: 10.1142/S0578563412500015)
- The special issue of Coastal Engineering Journal in 2012 contains many contributors of survey, please check individual papers for local use.


======================================================================
Data format
======================================================================

Data are text formatted with commna delimiter as follows. The asterisk "*" marked is important column of data set.
	#			serial number
	ID			ID for survey team
	location			measured location  (in Japanese)
	lon			longitude
	lat			latitude
	date			date of survey
	time			time of survey (JST = +9 UTC)
	measured height		* measured height 1 (raw data)
	height corrected by surveyor 	measured height 2 (tide corrected by individual survey team)
	height corrected by ttjt 	* measured height 3 (tide corrected by survey group) See Mori et al. (2012) CEJ 
	height from msl 		measured height 4 (height from MSL) See Mori et al. (2012) CEJ 
	height from TP 		measured height 5 (height from TP: Tokyo Peil) See Mori et al. (2012) CEJ 
	runup distance		distance from shoreline (unit:m)
	tide 1			tide level at survey  (individual survey team) (unit:m)
	tide 2			tide level at maximum tsunami arrived (individual survey team) (unit:m)
	tide 3			tide level at survey  (universally corrected by survey group) (unit:m)
	time 3			date and time for tide 3 (JST)
	tide 4			tide level at maximum tsunami arrived  (universally corrected by survey group) (unit:m)
	time 4			date and time for tide 4 (JST)
	type			type of data
	reliability			* reliability
	target			target of survey (in Japanese)
	mark			reason of survey (in Japanese)
	group			name of survey team  (in Japanese)
	comment		comment  (in Japanese)
	file type			method of  measurement
	flag for time		void
	original file name		original excel file name
- The tsunami height is generally defined elevation from astronomical tide when it arrived. 
- The unit of location is degree, the unit of height is meter and the unit of time is hour. The time is indicated fraction below the hour.
- The unused or not measured column is indicated "NaN" (Not a Number)
- There two different tsunami heights (inundation height) as
	1) Tide corrected by survey secretary entire dataset universally.
	2) Tide corrected by individual survey team
- The type 1) of tsunami heights are provided as relative height to astronomical height (measured height 3), Mean Sea Level (measured height 4) or Tokyo Peil (measured height 5).
- The type 2) of tsunami height is not shown for all of data.
- The inundation height is provided depends on survey team.
- The type of data are
	R: Run-up height
	I: inundation height
	P: mark in the harbor 
	W: very weak and could not measured
- Reliability
	A: High, clear mark, small error of measurement
	B: Middle, not clear mark but reliable information from witness, small error of measurement
	C: Low, run-up on the beach or location far from shoreline, large error of measurement
	D: Marginal, not clear, large uncertainty of measurement


======================================================================
Making survey dataset
======================================================================

1. Initial screening of data
2. Pre-process of tide correction
	The maximum tsunami arrival time was estimated by numerical simulations
	The astronomical tide near field was estimated by tidal data base by  National Astronomical Observatory of Japan
	The astronomical tide far field (west of Chiba and east of Hokkaido) was estimated by measured tide.
3. Tide was corrected based on item 2
4. The detail of making process has been summarized 
    Mori, N., T. Takahashi and The 2011 Tohoku Earthquake Tsunami Joint Survey Group (2012) Nationwide survey of the 2011 Tohoku earthquake tsunami, Coastal Engineering Journal, in press, Vol.54, Issue 1, pp.1-27. 


======================================================================
Known problems
======================================================================

- raw 537 and 538 are tentative
- duplicate points
	raw 3489 and 3497
	raw 3571 and 3580

======================================================================
 Contact
======================================================================

Any questions about the survey dataset and survey itself are welcome. 
- ML of survey team (restricted by subscriber)
	2011-ttjt@cm.kansai-u.ac.jp
- TTJT secretariat
	dpri@oceanwave.jp
We are a volunteer group, so please do not expect quick reply to your question.

======================================================================
 Update History
======================================================================

pre-release 20110520
	(only in Japanese)

release 20110621
	(only in Japanese)

release 20110715
	(only in Japanese)

release 20110826
	(only in Japanese)

release 20110907 (20110826c)
	(only in Japanese)

release 20110907 
	(only in Japanese)

release 20111017
	(only in Japanese)

release 20111031
	(only in Japanese)

release 20111110
	(only in Japanese)

release 20111201
	(only in Japanese)

release 20111215
	(only in Japanese)

release 20120101
	(only in Japanese)

release 20120203 Major
(We do not recommend to use earlier dataset before 20120203 due to  quality reason)
	The survey sheets by Miyagi Prefecture have been modified.

release 20120208 Major
	The data from Fukushima restricted area has been included.
	The data sheet by BRI has been modified.

release 20120301 Major
	The survey sheets by Miyagi Prefecture have been modified.
	The survey sheets by Nagaoka Institute of Technology have been modified.

release 20120330 Minor
	The locations of three data sheet by YNU have been modified.

release 201204 17 Minor
	The dataset has been slightly modified

release 20120330 Minor
	The locations of three data sheet by YNU have been modified.

release 20120425 Minor
	locations of 20 points were revised

release 20120806 Major
	39 points in Fukushima were added by UTokyo, Tokyo University of Marine Science and Technology, Geogia Tech. and Oregon State U.
	Two data set by NDA and DPRI-KU have been modified.

release 20120929 Major
	38 points in Fukushima were added by UTokyo, Tokyo University of Marine Science and Technology, Geogia Tech. and Oregon State U.

release 20121003
	6 points in Iwate by Geogia Tech have been added
	Total 5327 points

release 20121229
	Added data by Tokyo University of Science 
	Added data by National Institute for Land and Infrastructure Management
	Added data by PARI
	Added data by Tokai University
	Added data by  National Research Institute of Fisheries Engineering
	Total 5907 points
	This release will be final dataset from survey group.
