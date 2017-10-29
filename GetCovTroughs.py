#!/usr/bin/python

import numpy

__Njk__ = 30 

__quantile__ = [ '0', '1', '2', '3', '4' ]
__apperture__ = [ '10', '20', '30' ]

if __name__ == "__main__":
	wdd = {}
	cdd = {}
	wjk = {}
	wav = {}

	for band_ in ['mag_auto_i']:
		wdd[band_] = {}
		wav[band_] = {}
		cdd[band_] = {}
		wjk[band_] = {}

		for ap_ in __apperture__:
			for qt_ in __quantile__:
				for ii_ in xrange(5):
					wjk[band_][ii_] = {}
					cdd[band_][ii_] = numpy.zeros( (7,7) ) 
					wav[band_][ii_] = numpy.zeros( (2,7) )

					wdd[band_][ii_] = numpy.loadtxt('./results_number_counts_troughs/w_nc_troughs_ap'+str(ap_)+'_qt'+str(qt_)+'_'+band_+str(ii_),usecols=[0,3]).T
		
					for jk_ in xrange(__Njk__):
						wjk[band_][ii_][jk_] = numpy.loadtxt('./results_number_counts_troughs/w_nc_troughs_ap'+str(ap_)+'_qt'+str(qt_)+'_'+band_+str(ii_)+'_jk'+str(jk_),usecols=[0,3]).T

					for th_ in xrange(7):
						wav[band_][ii_][0][th_] = wjk[band_][ii_][0][0][th_]
						wav[band_][ii_][1][th_] = numpy.average( [ wjk[band_][ii_][jk_][1][th_] for jk_ in xrange(__Njk__) ] )

					for th1_ in xrange(7):
						for th2_ in xrange(7):
							for jk_ in xrange(__Njk__):
								cdd[band_][ii_][th1_][th2_] += (wdd[band_][ii_][1][th1_]-wjk[band_][ii_][jk_][1][th1_])*(wdd[band_][ii_][1][th2_]-wjk[band_][ii_][jk_][1][th2_])

					numpy.savetxt('./results_number_counts_troughs/cov_nc_ap'+ap_+'_qt'+qt_+'_'+band_+str(ii_),cdd[band_][ii_])
					numpy.savetxt('./results_number_counts_troughs/w_av_nc_ap'+ap_+'_qt'+qt_+'_'+band_+str(ii_),wdd[band_][ii_].T)
			
