#!/usr/bin/python

import numpy
import math
import treecorr
from astropy.io import fits
import matplotlib.pyplot as plt

__quantile__ = [ '0', '1', '2', '3', '4' ]
__apperture__ = [ '10', '20', '30', '60' ]

__bands__ = ['mag_auto_r', 'mag_auto_i', 'mag_auto_z']
__edges__ = { 'mag_auto_i':[ 20.0, 20.5, 21.0, 21.5, 22.0 ] }

__zbins__ = numpy.linspace(0,1.5,31)

if __name__ == "__main__":

	red = fits.open('y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10.fit')[1].data
	red = red[ (red['zspec'] > 0)  & (0.0 < red['zredmagic']) & (red['zredmagic'] < 0.45)]

	sour = numpy.hstack([fits.open('./processed_data/y1a1_gold_bpz_0.8_0.85.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_0.85_0.9.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_0.9_0.95.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_0.95_1.0.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.0_1.05.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.05_1.1.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.1_1.15.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.15_1.2.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.15_1.2.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.2_1.25.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.25_1.3.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.3_1.35.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.35_1.4.fits')[1].data, \
                             fits.open('./processed_data/y1a1_gold_bpz_1.45_1.5.fits')[1].data])

	sour = sour[ (sour['mlg'] > 0) & (sour['mlr'] > 0) & (sour['mli'] > 22.0) & (sour['mlz'] > 0) & (sour['flg'] >= 1)]
	sour = sour[ (-60 < sour['dec']) & (sour['dec'] < -40) ]

	sour = sour[(-1 < sour['mag_auto_g']-sour['mag_auto_r']) & (sour['mag_auto_g']-sour['mag_auto_r'] < 3) ]
	sour = sour[(-1 < sour['mag_auto_r']-sour['mag_auto_i']) & (sour['mag_auto_r']-sour['mag_auto_i'] < 2) ]
	sour = sour[(-1 < sour['mag_auto_i']-sour['mag_auto_z']) & (sour['mag_auto_i']-sour['mag_auto_z'] < 2) ]

	sour = sour[ sour['modest_class'] == 1]
	sour = sour[ sour['spread_model_i'] + (5./3.)*sour['spreaderr_model_i'] > 0.007 ]

	rans = fits.open('./processed_data/y1a1_random_uniform_masked.fits')[1].data

	rans = rans[ (rans['mlg'] > 0) & (rans['mlr'] > 0) & (rans['mli'] > 22.0) & (rans['mlz'] > 0) & (rans['flg'] >= 1)]
	rans = rans[ (-60 < rans['dec']) & (rans['dec'] < -40) ]
	ranl = rans

	lens = fits.open('./voids_masked.fits')[1].data

	lens = lens[ (-60 < lens['dec']) & (lens['dec'] < -40) ]
	lens = lens[ (lens['mlr'] > 0) & (lens['mli'] > 22.0) & (lens['mlz'] > 0) & (lens['good'] >= 1)]
	lens = lens[ lens['z'] < 0.45 ]


	for band_ in ['mag_auto_i']:
		for ii_ in xrange(5):

			print band_,ii_

			mySour = sour[ sour[band_] < __edges__[band_][ii_] ]

			plt.cla()
			plt.hist(mySour['z_mc_bpz'],bins=100,range=[0,2],normed=1,alpha=0.4)
			plt.hist(red['zspec'],bins=100,range=[0,2],normed=1,alpha=0.4)
			plt.savefig('./control_plots_voids/phiz_'+band_+str(ii_)+'.png')

			plt.cla()
			plt.hist2d(mySour['ra'],mySour['dec'],bins=[100,100])
			plt.savefig('./control_plots_voids/density_sour_'+band_+str(ii_)+'.png')

			lens_cat = treecorr.Catalog(ra=lens['ra'],dec=lens['dec'],ra_units='degrees',dec_units='degrees')
			ranl_cat = treecorr.Catalog(ra=ranl['ra'],dec=ranl['dec'],ra_units='degrees',dec_units='degrees')
			sour_cat = treecorr.Catalog(ra=mySour['ra'],dec=mySour['dec'],ra_units='degrees',dec_units='degrees')
			rans_cat = treecorr.Catalog(ra=rans['ra'],dec=rans['dec'],ra_units='degrees',dec_units='degrees')

			dd = treecorr.NNCorrelation(nbins=15,min_sep=0.1,max_sep=20,sep_units='degrees')
			dr = treecorr.NNCorrelation(nbins=15,min_sep=0.1,max_sep=20,sep_units='degrees')
			rd = treecorr.NNCorrelation(nbins=15,min_sep=0.1,max_sep=20,sep_units='degrees')
			rr = treecorr.NNCorrelation(nbins=15,min_sep=0.1,max_sep=20,sep_units='degrees')

			dd.process(lens_cat,sour_cat,num_threads=20)
			dr.process(lens_cat,rans_cat,num_threads=20)
			rd.process(ranl_cat,sour_cat,num_threads=20)
			rr.process(ranl_cat,rans_cat,num_threads=20)

			dd.write('./results_number_counts_voids/w_nc_troughs_'+band_+str(ii_),dr=dr,rd=rd,rr=rr)
