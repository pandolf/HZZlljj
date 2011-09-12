# Simple counting experiment, with one signal and one background process
#imax 1  number of channels
#jmax 1  number of backgrounds 
#kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
shapes ggH CMS_hzz2l2q_mm2b hzz2l2q_mm2b.input.root  w:signal
shapes VBF CMS_hzz2l2q_mm2b hzz2l2q_mm2b.input.root  w:signal 
shapes background CMS_hzz2l2q_mm2b hzz2l2q_mm2b.input.root w:background
###shapes data_obs CMS_hzz2l2q_mm2b hzz2l2q_mm2b.input.root w:data_obs
shapes data_obs CMS_hzz2l2q_mm2b hzz2l2q_mm2b.input.root w:dataset_obs
------------
bin         CMS_hzz2l2q_mm2b
observation <dummyobs>
------------
bin                CMS_hzz2l2q_mm2b	  CMS_hzz2l2q_mm2b	          CMS_hzz2l2q_mm2b
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
#CMS_recoe	lnN		1.015	1.015	1.0
CMS_scale_j	<dummyJES>
CMS_eff_b 	<dummybeff>
CMS_hzz2l2q_pu		lnN		1.02	        1.02		1.0			      	
#CMS_hzz2l2q_qgsep2b   lnN	 	<developing...>
#CMS_hzz2l2q_sig2bp0
#...                                                  <developing...>
#CMS_hzz2l2q_sig2bp5
CMS_hzz2l2q_bkg2bmmp0   <dummybnorm>
#CMS_hzz2l2q_bkg2bp1        param   226.23     ----
CMS_hzz2l2q_bkg2bp2         param  72.124    11.9277
#CMS_hzz2l2q_bkg2bp3      param    2.7583 ---
CMS_hzz2l2q_bkg2bp4         param  -0.457167   0.311

