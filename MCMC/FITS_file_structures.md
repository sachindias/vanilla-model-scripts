# FITS File Structure

For those who are able to access the data, it consists of several ~ 400MB files. Due to the complexities during warm-up and having to guide the walkers (to prevent them getting stuck in local minima) files have complicated names. These are detailed here for ease of understanding.

### Key:
| substring| meaning|
| ------ | ------ |
| afree | data where spin is a free parameter |
| a0 | data where spin is fixed at a=0|
| amax | data where spin fixed at a=0.998|


## GRO6:

rev0970_27_9_23_afree_10K.fits
rev0970_27_9_23_afree_2M.fits
rev0970_27_9_23_afree_2M_P2.fits
rev0970_27_9_23_afree_2M_P2_ext.fits
rev0970_27_9_23_afree_2M_P3.fits
rev0970_27_9_23_afree_2M_P3_ext.fits
rev0970_27_9_23_afree_2M_P3_ext2.fits
...
rev0970_27_9_23_afree_2M_P3_ext42.fits
rev0970_27_9_23_afree_2M_P3_ext42_P2.fits
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext.fits
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext2.fits
...
rev0970_27_9_23_afree_2M_P3_ext42_P2_ext103.fits