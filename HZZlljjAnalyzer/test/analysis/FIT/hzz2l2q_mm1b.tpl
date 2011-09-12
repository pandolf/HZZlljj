# Simple counting experiment, with one signal and one background process
#imax 1  number of channels
#jmax 1  number of backgrounds 
#kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
shapes ggH CMS_hzz2l2q_mm1b hzz2l2q_mm1b.input.root  w:signal
shapes VBF CMS_hzz2l2q_mm1b hzz2l2q_mm1b.input.root  w:signal 
shapes background CMS_hzz2l2q_mm1b hzz2l2q_mm1b.input.root w:background
###shapes data_obs CMS_hzz2l2q_mm1b hzz2l2q_mm1b.input.root w:data_obs
shapes data_obs CMS_hzz2l2q_mm1b hzz2l2q_mm1b.input.root w:dataset_obs
------------
# we have just one channel, in which we observe 0 events
bin         CMS_hzz2l2q_mm1b
observation <dummyobs>
------------
bin                CMS_hzz2l2q_mm1b	  CMS_hzz2l2q_mm1b	          CMS_hzz2l2q_mm1b
process       ggH       		  				  VBF                 background
process         -1                                                          0                        1        
rate	    <dummy1>
------------
lumi		lnN	1.045			1.045			1.0
pdf_ggH	         <dummypdfggH>
pdf_qqH	       <dummypdfqqH>
QCDscale_ggH				 <dummyggH>
QCDscale_qqH			 <dummyVBF>
#theory_gamma                      <dummygammaBW>
CMS_trigger_m	lnN	1.02	1.02	1.0	
CMS_eff_m	lnN	1.008	1.008	1.0
CMS_scale_m	lnN	1.01	1.01	1.0
#CMS_recom	lnN	1.015	1.015	1.0
CMS_scale_j	<dummyJES>
CMS_eff_b 	<dummybeff>
CMS_hzz2l2q_pu		lnN		1.02	        1.02		1.0			      	
#CMS_hzz2l2q_qgsep1b   lnN	 	<developing...>
#CMS_hzz2l2q_sig1bp0
#...                                                  <developing...>
#CMS_hzz2l2q_sig1bp5
CMS_hzz2l2q_bkg1bmmp0   <dummybnorm>
#CMS_hzz2l2q_bkg1bp1        param  166.6  ---
CMS_hzz2l2q_bkg1bp2         param  87.4666 7.03731
#CMS_hzz2l2q_bkg1bp3      param   21.499 ---  
CMS_hzz2l2q_bkg1bp4         param  0.249461  0.0552891
