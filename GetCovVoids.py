#!/usr/bin/python

import numpy

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

		for ii_ in xrange(5):
			wjk[band_][ii_] = {}
			cdd[band_][ii_] = numpy.zeros( (15,15) ) 
			wav[band_][ii_] = numpy.zeros( (2,15) )

			wdd[band_][ii_] = numpy.loadtxt('./results_number_counts_voids/w_nc_troughs_'+band_+str(ii_),usecols=[0,3]).T

			for jk_ in xrange(192):
				wjk[band_][ii_][jk_] = numpy.loadtxt('./results_number_counts_voids/w_nc_'+band_+str(ii_)+'_jk'+str(jk_),usecols=[0,3]).T

			for th_ in xrange(15):
				wav[band_][ii_][0][th_] = wjk[band_][ii_][0][0][th_]
				wav[band_][ii_][1][th_] = numpy.average( [ wjk[band_][ii_][jk_][1][th_] for jk_ in xrange(192) ] )

			for th1_ in xrange(15):
				for th2_ in xrange(15):
					for jk_ in xrange(192):
						cdd[band_][ii_][th1_][th2_] += (wdd[band_][ii_][1][th1_]-wjk[band_][ii_][jk_][1][th1_])*(wdd[band_][ii_][1][th2_]-wjk[band_][ii_][jk_][1][th2_])

			numpy.savetxt('./results_number_counts_voids/cov_nc_'+band_+str(ii_),cdd[band_][ii_])
			numpy.savetxt('./results_number_counts_voids/w_av_nc_'+band_+str(ii_),wdd[band_][ii_].T)
			
